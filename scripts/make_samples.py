#!/usr/bin/env python3
"""Build samples.tsv from the CDC SARS-CoV-2 benchmark dataset TSVs.

Usage: python make_samples.py sars-cov-2-voivoc.tsv sars-cov-2-nonvoivoc.tsv sars-cov-2-failedQC.tsv
"""
import csv, sys

GROUPS = {"voivoc": "voivoc", "nonvoivoc": "nonvoivoc", "failedqc": "failedqc"}

def group_for(path):
    name = path.lower().rsplit("/", 1)[-1].replace("sars-cov-2-", "").replace(".tsv", "")
    return GROUPS[name]

rows = []
for path in sys.argv[1:]:
    group = group_for(path)
    with open(path, newline="") as fh:
        lines = list(csv.reader(fh, delimiter="\t"))
    # The file has a metadata header block; the real table starts at the row containing SRArun_acc
    start = next(i for i, l in enumerate(lines) if "SRArun_acc" in l)
    header = lines[start]
    acc_i = header.index("SRArun_acc")
    lin_i = next(i for i, h in enumerate(header) if h.lower().startswith("lineage"))
    note_i = header.index("benchmark_dataset") if "benchmark_dataset" in header else None
    for l in lines[start + 1:]:
        if len(l) <= acc_i or not l[acc_i] or l[acc_i] in ("NA", "-"):
            continue  # skips the Wuhan-1 reference row and blank lines
        lineage = l[lin_i] if l[lin_i] not in ("", "-") else "NA"
        note = l[note_i] if note_i is not None else ""
        rows.append((l[acc_i], group, lineage, note))

with open("samples.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "group", "lineage", "note"])
    w.writerows(rows)

from collections import Counter
print(f"Wrote samples.tsv with {len(rows)} samples:", dict(Counter(r[1] for r in rows)))