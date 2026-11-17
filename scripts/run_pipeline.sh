#!/usr/bin/env bash
# Reads-to-variants pipeline for all samples in samples.tsv
#   fastp -> bwa mem -> samtools sort/index -> QC -> bcftools mpileup/call -> filter
#
# Run from the project folder:   bash scripts/run_pipeline.sh
# Re-running is safe: samples whose final VCF already exists are skipped.

set -uo pipefail

REF=data/ref/NC_045512.2.fasta
THREADS=4
MIN_QUAL=30        # variant quality
MIN_DP=20          # reads covering the site
MIN_AF=0.75        # fraction of reads supporting the variant (single infection ~1.0)
MIN_GAP_DEPTH=10   # positions below this depth count as "uncovered" for the gap report

mkdir -p data/trim data/bam results/qc results/vcf

# warn if disk space is low (each sample needs up to ~1 GB of temporary space while it runs)
free_gb=$(df -g . | awk 'NR==2 {print $4}')
echo "Free disk space: ${free_gb} GB"
if [ "${free_gb:-0}" -lt 5 ]; then echo "Less than 5 GB free - free up space before running."; exit 1; fi
SUMMARY=results/qc/summary.tsv
[ -f "$SUMMARY" ] || printf "accession\tgroup\tlineage\tread_pairs\tpct_mapped\tpct_properly_paired\tpct_genome_covered\tmean_depth\tlongest_gap_bp\tn_variants\tstatus\n" > "$SUMMARY"

tail -n +2 samples.tsv | while IFS=$'\t' read -r acc group lineage note; do
  if [ -f results/vcf/${acc}.filtered.vcf.gz ]; then echo "[$acc] already done, skipping"; continue; fi
  echo "[$acc] starting ($group, $lineage)"

  # --- find the raw reads (.fastq or .fastq.gz; paired or single-end) ---
  r1=$(ls data/raw/${acc}_1.fastq* 2>/dev/null | head -1)
  r2=$(ls data/raw/${acc}_2.fastq* 2>/dev/null | head -1)
  [ -z "$r1" ] && r1=$(ls data/raw/${acc}.fastq* 2>/dev/null | head -1)
  if [ -z "$r1" ]; then
    echo "[$acc] no reads found, skipping"
    printf "%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tno_reads\n" "$acc" "$group" "$lineage" >> "$SUMMARY"
    continue
  fi

  # --- 1. trim + 2. align ---
  if [ -n "$r2" ]; then
    fastp -i "$r1" -I "$r2" -o data/trim/${acc}_1.fq.gz -O data/trim/${acc}_2.fq.gz \
          -j results/qc/${acc}.fastp.json -h results/qc/${acc}.fastp.html -w $THREADS 2> results/qc/${acc}.fastp.log
    bwa mem -t $THREADS -R "@RG\tID:${acc}\tSM:${acc}" $REF data/trim/${acc}_1.fq.gz data/trim/${acc}_2.fq.gz 2> results/qc/${acc}.bwa.log \
      | samtools sort -o data/bam/${acc}.bam -
  else
    echo "[$acc] single-end reads"
    fastp -i "$r1" -o data/trim/${acc}.fq.gz \
          -j results/qc/${acc}.fastp.json -h results/qc/${acc}.fastp.html -w $THREADS 2> results/qc/${acc}.fastp.log
    bwa mem -t $THREADS -R "@RG\tID:${acc}\tSM:${acc}" $REF data/trim/${acc}.fq.gz 2> results/qc/${acc}.bwa.log \
      | samtools sort -o data/bam/${acc}.bam -
  fi
  # stop here if alignment failed (e.g. disk full), and record it instead of reporting empty results as "ok"
  if [ ! -s data/bam/${acc}.bam ] || ! samtools quickcheck data/bam/${acc}.bam; then
    echo "[$acc] ALIGNMENT FAILED, see results/qc/${acc}.bwa.log"
    rm -f data/bam/${acc}.bam data/trim/${acc}*
    printf "%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\talign_failed\n" "$acc" "$group" "$lineage" >> "$SUMMARY"
    continue
  fi
  samtools index data/bam/${acc}.bam
  rm -f data/trim/${acc}*   # trimmed reads are no longer needed once the BAM exists; saves several GB

  # --- 3. alignment QC ---
  samtools flagstat data/bam/${acc}.bam > results/qc/${acc}.flagstat.txt
  samtools coverage data/bam/${acc}.bam > results/qc/${acc}.coverage.txt
  pairs=$(grep -m1 " read1" results/qc/${acc}.flagstat.txt | cut -d' ' -f1)
  pct_mapped=$(grep -m1 "primary mapped (" results/qc/${acc}.flagstat.txt | sed -E 's/.*\(([0-9.]+)%.*/\1/')
  pct_pp=$(grep -m1 "properly paired (" results/qc/${acc}.flagstat.txt | sed -E 's/.*\(([0-9.]+)%.*/\1/')
  [ -z "$pct_pp" ] && pct_pp=NA
  read -r pct_cov mean_depth < <(awk 'NR==2 {print $6, $7}' results/qc/${acc}.coverage.txt)
  longest_gap=$(samtools depth -a data/bam/${acc}.bam \
    | awk -v m=$MIN_GAP_DEPTH '$3<m {if(s && $2==e+1){e=$2} else {if(s && e-s+1>best) best=e-s+1; s=$2; e=$2}} END{if(s && e-s+1>best) best=e-s+1; print best+0}')

  # --- 4. call variants ---
  #   -L 10000 : raise indel-calling depth limit (default 250 is too low for amplicon data)
  bcftools mpileup -f $REF -d 10000 -L 10000 -a AD,DP -Ou data/bam/${acc}.bam 2> results/qc/${acc}.mpileup.log \
    | bcftools call --ploidy 1 -mv -Oz -o results/vcf/${acc}.raw.vcf.gz
  bcftools index -f results/vcf/${acc}.raw.vcf.gz

  # --- 5. indel check against ALL reads ---
  #   For indels, bcftools' AD/DP count only the reads it re-aligned around the candidate, which can make a
  #   rare artefact look common (e.g. 123/142 reads instead of 121/6963). So for every indel we count, from the
  #   raw pileup, how many of ALL reads at the site carry exactly that indel, and exclude it if below MIN_AF.
  : > results/vcf/${acc}.indel_exclude.tsv
  : > results/qc/${acc}.indel_af.tsv
  bcftools query -i 'TYPE="indel"' -f '%CHROM\t%POS\t%REF\t%ALT\n' results/vcf/${acc}.raw.vcf.gz \
  | while IFS=$'\t' read -r chrom pos ref alt; do
      len=$(( ${#alt} - ${#ref} ))                       # negative = deletion, positive = insertion
      if [ $len -lt 0 ]; then sign="-"; n=$(( -len )); else sign="+"; n=$len; fi
      samtools mpileup -r "${chrom}:${pos}-${pos}" -f $REF -d 0 -Q 0 data/bam/${acc}.bam 2>/dev/null \
      | awk -v s="$sign" -v n="$n" -v chrom="$chrom" -v pos="$pos" -v min="$MIN_AF" \
            -v ex="results/vcf/${acc}.indel_exclude.tsv" -v logf="results/qc/${acc}.indel_af.tsv" '
          { pat = "\\" s n "[ACGTNacgtn]"; k = gsub(pat, "", $5); af = ($4 > 0) ? k / $4 : 0
            printf "%s\t%s\t%s%s\t%d\t%d\t%.4f\t%s\n", chrom, pos, s, n, k, $4, af, (af >= min ? "keep" : "exclude") >> logf
            if (af < min) print chrom "\t" pos >> ex }'
    done

  # --- 6. filter: quality, depth, allele frequency (SNPs via AD), and the indel check above ---
  exclude_opt=""
  [ -s results/vcf/${acc}.indel_exclude.tsv ] && exclude_opt="-T ^results/vcf/${acc}.indel_exclude.tsv"
  bcftools view $exclude_opt -i "QUAL>=${MIN_QUAL} && INFO/DP>=${MIN_DP} && FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1])>=${MIN_AF}" \
    results/vcf/${acc}.raw.vcf.gz -Oz -o results/vcf/${acc}.filtered.vcf.gz
  bcftools index -f results/vcf/${acc}.filtered.vcf.gz
  if [ ! -s results/vcf/${acc}.filtered.vcf.gz ]; then
    echo "[$acc] VARIANT CALLING FAILED, see results/qc/${acc}.mpileup.log"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tvcf_failed\n" \
      "$acc" "$group" "$lineage" "$pairs" "$pct_mapped" "$pct_pp" "$pct_cov" "$mean_depth" "$longest_gap" >> "$SUMMARY"
    continue
  fi
  n_var=$(bcftools view -H results/vcf/${acc}.filtered.vcf.gz | wc -l | tr -d ' ')

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tok\n" \
    "$acc" "$group" "$lineage" "$pairs" "$pct_mapped" "$pct_pp" "$pct_cov" "$mean_depth" "$longest_gap" "$n_var" >> "$SUMMARY"
  echo "[$acc] done: ${pct_cov}% covered, mean depth ${mean_depth}, longest gap ${longest_gap} bp, ${n_var} variants"
done

echo "All samples processed. QC summary: $SUMMARY"
