#!/usr/bin/env python3
"""Build a samples x variants matrix from the filtered VCFs.

Run from the project folder:   python scripts/build_matrix.py

Encoding for each sample and variant:
    1   the variant was called in this sample
    0   not called, AND the site had enough reads (>= MIN_DEPTH) to have called it -> confidently reference
    NA  not called, but depth at the site was too low to tell (e.g. amplicon dropout)

Treating low-coverage sites as missing (not as 0) matters: otherwise a dropout looks like
"this sample lacks the mutation", which is exactly the trap seen in ERR5556343 (ORF8/N mutations).

Outputs (results/matrix/):
    variant_matrix.tsv    rows = samples, columns = variants (POS + REF > ALT)
    sample_metadata.tsv   group, lineage, QC values and whether the sample passes QC
    variant_summary.tsv   how often each variant occurs in each group
"""
import csv, gzip, os, subprocess, sys
from collections import defaultdict

MIN_DEPTH = 20          # reads needed at a site to call it confidently "reference"
MIN_COVERAGE = 90.0     # % genome covered for a sample to pass QC
GENOME_LEN = 29903

VCF_DIR, BAM_DIR = "results/vcf", "data/bam"
DEPTH_DIR, OUT_DIR = "results/depth", "results/matrix"
os.makedirs(DEPTH_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def depth_profile(acc):
    """Per-position depth (index 1..29903). Cached in results/depth so samtools only runs once per sample."""
    path = os.path.join(DEPTH_DIR, f"{acc}.depth.tsv.gz")
    if not os.path.exists(path):
        bam = os.path.join(BAM_DIR, f"{acc}.bam")
        out = subprocess.run(["samtools", "depth", "-a", bam], capture_output=True, text=True, check=True).stdout
        with gzip.open(path, "wt") as fh:
            fh.write(out)
    depth = [0] * (GENOME_LEN + 1)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            _, pos, d = line.split("\t")
            pos = int(pos)
            if pos <= GENOME_LEN:
                depth[pos] = int(d)
    return depth


def read_variants(acc):
    """Set of variant keys like '23063A>T' from the filtered VCF."""
    found = set()
    with gzip.open(os.path.join(VCF_DIR, f"{acc}.filtered.vcf.gz"), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            pos, ref = f[1], f[3]
            for alt in f[4].split(","):
                found.add((int(pos), ref, alt))
    return found


# ---- 1. samples and QC --------------------------------------------------------------
samples = {r["accession"]: r for r in read_tsv("samples.tsv")}
summary = {r["accession"]: r for r in read_tsv("results/qc/summary.tsv")}

usable = []
for acc, s in samples.items():
    q = summary.get(acc)
    if not q or q["status"] != "ok" or not os.path.exists(os.path.join(VCF_DIR, f"{acc}.filtered.vcf.gz")):
        print(f"  skipping {acc}: no finished pipeline output")
        continue
    usable.append(acc)

# ---- 2. collect variants ------------------------------------------------------------
calls = {acc: read_variants(acc) for acc in usable}
all_variants = sorted(set().union(*calls.values()))
print(f"{len(usable)} samples, {len(all_variants)} distinct variants")

# ---- 3. fill the matrix (1 / 0 / NA) ------------------------------------------------
matrix = {}
for i, acc in enumerate(usable, 1):
    depth = depth_profile(acc)
    row = []
    for v in all_variants:
        pos = v[0]
        if v in calls[acc]:
            row.append("1")
        elif depth[pos] >= MIN_DEPTH:
            row.append("0")
        else:
            row.append("NA")
    matrix[acc] = row
    print(f"  [{i}/{len(usable)}] {acc}: {row.count('1')} variants, {row.count('NA')} unknown sites")

labels = [f"{p}{r}>{a}" for p, r, a in all_variants]

with open(os.path.join(OUT_DIR, "variant_matrix.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["accession"] + labels)
    for acc in usable:
        w.writerow([acc] + matrix[acc])

# ---- 4. sample metadata (incl. QC pass/fail) ----------------------------------------
with open(os.path.join(OUT_DIR, "sample_metadata.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "group", "lineage", "pct_genome_covered", "mean_depth", "longest_gap_bp",
                "n_variants", "n_unknown_sites", "passes_qc"])
    for acc in usable:
        s, q = samples[acc], summary[acc]
        cov = float(q["pct_genome_covered"])
        w.writerow([acc, s["group"], s["lineage"], q["pct_genome_covered"], q["mean_depth"], q["longest_gap_bp"],
                    matrix[acc].count("1"), matrix[acc].count("NA"), "yes" if cov >= MIN_COVERAGE else "no"])

# ---- 5. per-variant summary by group ------------------------------------------------
groups = sorted({samples[a]["group"] for a in usable})
with open(os.path.join(OUT_DIR, "variant_summary.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    header = ["variant", "position", "n_samples_with", "n_unknown"]
    for g in groups:
        header += [f"{g}_with", f"{g}_known", f"{g}_freq"]
    w.writerow(header)
    for j, label in enumerate(labels):
        col = {acc: matrix[acc][j] for acc in usable}
        row = [label, all_variants[j][0], sum(v == "1" for v in col.values()), sum(v == "NA" for v in col.values())]
        for g in groups:
            vals = [col[a] for a in usable if samples[a]["group"] == g]
            known = [v for v in vals if v != "NA"]
            n_with = known.count("1")
            row += [n_with, len(known), f"{n_with / len(known):.2f}" if known else "NA"]
        w.writerow(row)

print(f"\nWrote {OUT_DIR}/variant_matrix.tsv, sample_metadata.tsv, variant_summary.tsv")
