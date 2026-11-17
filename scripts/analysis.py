#!/usr/bin/env python3
"""Can variant profiles distinguish VOI/VOC samples from other lineages?

Run from the project folder:   python scripts/analysis.py
Reads results/matrix/*.tsv, writes figures and tables to results/analysis/.

Steps
  1. Select samples: VOI/VOC + non-VOI/VOC that pass QC (>= 90% genome covered).
  2. Select features: variants seen in >= 2 of those samples (singletons can't generalise).
  3. Unsupervised: PCA + K-means (k=2) - do samples separate without labels?
  4. Supervised: Random Forest, evaluated two ways
       a) stratified 5-fold CV (repeated)       -> optimistic: same lineage can be in train and test
       b) leave-lineage-out (grouped CV)        -> honest: can it recognise a lineage it has never seen?
  5. Interpretation: top variants by permutation importance, mapped to genes.
  6. Robustness: retrain without the top 10 variants; predict the 2 VOI/VOC samples that failed QC.

Missing data (NA = site not covered) is filled with "0" INSIDE each training fold (SimpleImputer in the
pipeline), so no information from test samples leaks into training. The share of NA per sample is reported.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import (balanced_accuracy_score, confusion_matrix, adjusted_rand_score,
                             recall_score)
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedGroupKFold, cross_val_predict
from sklearn.inspection import permutation_importance

SEED = 42
MIN_SAMPLES_PER_VARIANT = 2
OUT = "results/analysis"
os.makedirs(OUT, exist_ok=True)
report = []
def log(msg=""):
    print(msg); report.append(msg)

# SARS-CoV-2 genes on NC_045512.2 (start, end)
GENES = [("5'UTR", 1, 265), ("ORF1ab", 266, 21555), ("S", 21563, 25384), ("ORF3a", 25393, 26220),
         ("E", 26245, 26472), ("M", 26523, 27191), ("ORF6", 27202, 27387), ("ORF7a", 27394, 27759),
         ("ORF7b", 27756, 27887), ("ORF8", 27894, 28259), ("N", 28274, 29533), ("ORF10", 29558, 29674),
         ("3'UTR", 29675, 29903)]
def gene_of(pos):
    hits = [g for g, s, e in GENES if s <= pos <= e]
    return "/".join(hits) if hits else "intergenic"

# ---- 1. data ----------------------------------------------------------------------------------
X_all = pd.read_csv("results/matrix/variant_matrix.tsv", sep="\t", index_col=0, na_values="NA")
meta = pd.read_csv("results/matrix/sample_metadata.tsv", sep="\t", index_col=0).loc[X_all.index]

keep = meta["group"].isin(["voivoc", "nonvoivoc"]) & (meta["passes_qc"] == "yes")
X, m = X_all[keep], meta[keep]
y = (m["group"] == "voivoc").astype(int).values
groups = m["lineage"].values

log("=== Samples ===")
log(f"Included: {len(m)} ({y.sum()} VOI/VOC, {(y == 0).sum()} non-VOI/VOC); "
    f"{m['lineage'].nunique()} lineages")
excluded = meta[meta["group"].isin(["voivoc", "nonvoivoc"]) & (meta["passes_qc"] != "yes")]
log(f"Excluded for QC (<90% covered): {', '.join(f'{a} ({r.lineage}, {r.pct_genome_covered:.0f}%)' for a, r in excluded.iterrows())}")

# ---- 2. features ------------------------------------------------------------------------------
counts = (X == 1).sum()
X = X.loc[:, counts >= MIN_SAMPLES_PER_VARIANT]
log(f"\n=== Features ===")
log(f"{X_all.shape[1]} variants in total; {X.shape[1]} seen in >= {MIN_SAMPLES_PER_VARIANT} included samples")
log(f"Missing (NA) cells among included samples: {X.isna().values.mean():.1%}")

# ---- 3. unsupervised --------------------------------------------------------------------------
Xi = SimpleImputer(strategy="constant", fill_value=0).fit_transform(X)
pca = PCA(n_components=2, random_state=SEED)
pcs = pca.fit_transform(Xi)
km = KMeans(n_clusters=2, n_init=20, random_state=SEED).fit(Xi)
log(f"\n=== Unsupervised ===")
log(f"PCA variance explained: PC1 {pca.explained_variance_ratio_[0]:.1%}, PC2 {pca.explained_variance_ratio_[1]:.1%}")
log(f"K-means (k=2) vs true groups: adjusted Rand index = {adjusted_rand_score(y, km.labels_):.2f} "
    f"(1 = perfect match, 0 = no better than chance)")

fig, ax = plt.subplots(figsize=(7, 5.5))
for label, colour, name in [(0, "#8a8a8a", "non-VOI/VOC"), (1, "#c0392b", "VOI/VOC")]:
    sel = y == label
    ax.scatter(pcs[sel, 0], pcs[sel, 1], c=colour, label=name, s=45, alpha=0.8, edgecolor="white")
for i, (acc, row) in enumerate(m.iterrows()):
    if y[i] == 1:
        ax.annotate(row["lineage"], pcs[i], fontsize=7, xytext=(4, 3), textcoords="offset points")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.0%} of variance)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.0%} of variance)")
ax.set_title("Samples by variant profile (PCA)")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/pca.png", dpi=200); plt.close(fig)

# ---- 4. supervised ----------------------------------------------------------------------------
def model():
    return make_pipeline(SimpleImputer(strategy="constant", fill_value=0),
                         RandomForestClassifier(n_estimators=500, class_weight="balanced", random_state=SEED))

def evaluate(Xd, label):
    """Returns (repeated-stratified scores, leave-lineage-out predictions)."""
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=SEED)
    bal, rec_voc = [], []
    for tr, te in rskf.split(Xd, y):
        clf = model().fit(Xd.iloc[tr], y[tr]); p = clf.predict(Xd.iloc[te])
        bal.append(balanced_accuracy_score(y[te], p)); rec_voc.append(recall_score(y[te], p))
    sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED)
    p_lineage = cross_val_predict(model(), Xd, y, cv=sgkf, groups=groups)
    log(f"{label}")
    log(f"  stratified 5-fold x10   balanced accuracy {np.mean(bal):.2f} +/- {np.std(bal):.2f} | "
        f"VOI/VOC recall {np.mean(rec_voc):.2f}")
    log(f"  leave-lineage-out       balanced accuracy {balanced_accuracy_score(y, p_lineage):.2f} | "
        f"VOI/VOC recall {recall_score(y, p_lineage):.2f}")
    return p_lineage

log(f"\n=== Classification (VOI/VOC vs non-VOI/VOC) ===")
dummy = cross_val_predict(DummyClassifier(strategy="stratified", random_state=SEED), X.fillna(0), y,
                          cv=StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED), groups=groups)
log(f"Baseline (random guessing in class proportions): balanced accuracy {balanced_accuracy_score(y, dummy):.2f}")
p_lin = evaluate(X, "Random Forest, all variants:")

cm = confusion_matrix(y, p_lin)
fig, ax = plt.subplots(figsize=(4, 3.6))
ax.imshow(cm, cmap="Blues")
for (i, j), v in np.ndenumerate(cm):
    ax.text(j, i, v, ha="center", va="center", color="white" if v > cm.max() / 2 else "black", fontsize=13)
ax.set_xticks([0, 1], ["non-VOI/VOC", "VOI/VOC"]); ax.set_yticks([0, 1], ["non-VOI/VOC", "VOI/VOC"])
ax.set_xlabel("Predicted"); ax.set_ylabel("True"); ax.set_title("Leave-lineage-out predictions")
fig.tight_layout(); fig.savefig(f"{OUT}/confusion_matrix_leave_lineage_out.png", dpi=200); plt.close(fig)

pred = m[["lineage"]].copy(); pred["true"] = np.where(y == 1, "VOI/VOC", "non-VOI/VOC")
pred["predicted_leave_lineage_out"] = np.where(p_lin == 1, "VOI/VOC", "non-VOI/VOC")
pred.to_csv(f"{OUT}/predictions_leave_lineage_out.tsv", sep="\t")
wrong = pred[pred["true"] != pred["predicted_leave_lineage_out"]]
log(f"  misclassified when lineage unseen: " + (", ".join(f"{a} ({r.lineage}, true {r.true})" for a, r in wrong.iterrows()) or "none"))

# ---- 5. interpretation ------------------------------------------------------------------------
final = model().fit(X, y)
imp = permutation_importance(final, X, y, n_repeats=30, random_state=SEED, scoring="balanced_accuracy")
rf_imp = final[-1].feature_importances_
tab = pd.DataFrame({"variant": X.columns, "position": [int("".join(c for c in v.split(">")[0] if c.isdigit())) for v in X.columns],
                    "rf_importance": rf_imp, "permutation_importance": imp.importances_mean})
tab["gene"] = tab["position"].map(gene_of)
tab["voivoc_freq"] = [(X.loc[y == 1, v] == 1).sum() / X.loc[y == 1, v].notna().sum() for v in X.columns]
tab["nonvoivoc_freq"] = [(X.loc[y == 0, v] == 1).sum() / max(X.loc[y == 0, v].notna().sum(), 1) for v in X.columns]
tab = tab.sort_values("rf_importance", ascending=False)
tab.to_csv(f"{OUT}/feature_importance.tsv", sep="\t", index=False, float_format="%.4f")
log(f"\n=== Top 15 variants (Random Forest importance) ===")
log(f"{'variant':24s}{'gene':8s}{'RF imp':>8s}{'VOC freq':>10s}{'non freq':>10s}")
for _, r in tab.head(15).iterrows():
    log(f"{r.variant:24s}{r.gene:8s}{r.rf_importance:8.3f}{r.voivoc_freq:10.2f}{r.nonvoivoc_freq:10.2f}")

top = tab.head(15)[::-1]
fig, ax = plt.subplots(figsize=(6.5, 5))
ax.barh([f"{v} ({g})" for v, g in zip(top.variant, top.gene)], top.rf_importance, color="#c0392b")
ax.set_xlabel("Random Forest importance"); ax.set_title("Variants that best separate VOI/VOC")
fig.tight_layout(); fig.savefig(f"{OUT}/feature_importance.png", dpi=200); plt.close(fig)

# ---- 6. robustness ----------------------------------------------------------------------------
log(f"\n=== Robustness ===")
X_drop = X.drop(columns=tab.head(10)["variant"])
evaluate(X_drop, "Random Forest, top 10 variants removed:")

if len(excluded):
    X_ex = X_all.loc[excluded.index, X.columns]
    p_ex = final.predict(X_ex)
    log("VOI/VOC samples that failed QC, predicted by the final model:")
    for (acc, r), p, na in zip(excluded.iterrows(), p_ex, X_ex.isna().mean(axis=1)):
        log(f"  {acc} ({r.lineage}, {r.pct_genome_covered:.0f}% covered, {na:.0%} of features unknown) "
            f"-> {'VOI/VOC' if p == 1 else 'non-VOI/VOC'}")

with open(f"{OUT}/results.txt", "w") as fh:
    fh.write("\n".join(report) + "\n")
print(f"\nSaved figures and tables to {OUT}/")
