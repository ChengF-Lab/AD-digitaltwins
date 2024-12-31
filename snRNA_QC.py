### Load package and read data

# python3.10
import scanpy as sc
import logging
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import scipy as sp
from scipy import sparse
from scipy.stats import median_abs_deviation
import doubletdetection
import scrublet as scr


adata = sc.read_h5ad("SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
adata


### 1. check QC


## mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
## ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
## hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
## QC metric using calculate_qc_metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=None, log1p=False)
adata

SHAPES = [["RAW", *adata.shape], ]


### 2. Plot raw quality
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], save="_rawQC_counts_plot.pdf")  

sns.histplot(adata.obs["pct_counts_mt"])
plt.savefig("rawQC1_plot.pdf")
plt.close()

sns.displot(adata.obs["total_counts"], bins=100, kde=False)
plt.savefig("rawQC2_displot_total_counts.pdf")
plt.close()

sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=True,save="_rawQC3_plot.pdf")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=True,save="_rawQC4_plot.pdf")

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",save="_rawQC5_plot.pdf")

sns.histplot(adata.obs["pct_counts_ribo"])
plt.savefig("rawQC6_plot.pdf")
plt.close()




### 3. remove reference

adata = adata[adata.obs["Cognitive Status"] != "Reference"]
SHAPES.append(["DATASET_FILTERING", *adata.shape])
SHAPES


### 4. Universal QC

# Universal QC
QCs = [
    adata.var["n_cells_by_counts"] >= 10,
    adata.obs["n_genes_by_counts"] >= 200,
    adata.obs["total_counts"] > 500,
]
QC_ROW_PASS = np.ones(adata.shape[0], dtype=bool) & QCs[1] & QCs[2]
QC_COL_PASS = np.ones(adata.shape[1], dtype=bool) & QCs[0]

SHAPES.append(["QC_1_PASS", QC_ROW_PASS.sum(), QC_COL_PASS.sum()])
SHAPES.append(["  MIN_CELLS>=10", "-", QCs[0].sum()])
SHAPES.append(["  MIN_GENES>=200", QCs[1].sum(), "-"])
SHAPES.append(["  TOTAL_COUNTS>500", QCs[2].sum(), "-"])
SHAPES
adata = adata[QC_ROW_PASS, :][:, QC_COL_PASS]



### 6. Adaptive QC

QCs = []
for field in ("total_counts", "n_genes_by_counts"):
    values = adata.obs[field]
    q75, q25 = np.percentile(values, [75, 25])
    lower_co = q25 - 3 * (q75 - q25)
    upper_co = q75 + 3 * (q75 - q25)
    QCs.append((values > lower_co) & (values < upper_co))
    print(f"Threshold for {field}: {lower_co} < x < {upper_co}")

for field, threshold in (
        ("pct_counts_mt", 10),
        ("pct_counts_rb", 20)
):
    values = adata.obs[field]
    q75, q25 = np.percentile(values, [75, 25])
    upper_co = q75 + 3 * (q75 - q25)
    QCs.append(values < max(upper_co, threshold))
    print(f"Threshold for {field}: x < {max(upper_co, threshold)}")

QC_PASS = np.ones(adata.shape[0], dtype=bool) & QCs[0] & QCs[1] & QCs[2] & QCs[3]

SHAPES.append(["QC_2_PASS", QC_PASS.sum(), "-"])
SHAPES.append(["  TOTAL_COUNTS<3IQR", QCs[0].sum(), "-"])
SHAPES.append(["  N_GENES_BY_COUNTS<3IQR", QCs[1].sum(), "-"])
SHAPES.append(["  MT_PERCENT<3IQR", QCs[2].sum(), "-"])
SHAPES.append(["  RB_PERCENT<3IQR", QCs[3].sum(), "-"])

adata = adata[QC_PASS]

SHAPES



### 7. Doublet detection

datasets = []
for sample in adata.obs["Donor ID"].unique():
    sample_adata = adata[adata.obs["Donor ID"] == sample].copy()
    expected = sample_adata.shape[0] / 1000 * 0.008
    print("Expected doublet rate for sample:", sample, expected)
    scrub = scr.Scrublet(sample_adata.X, expected_doublet_rate=expected)
    sample_adata.obs["doublet_scores"], sample_adata.obs["predicted_doublets"] = scrub.scrub_doublets(n_prin_comps=30)
    datasets.append(sample_adata)

# Doublet removal
adata = sc.concat(datasets, merge="unique")
adata = adata[adata.obs["predicted_doublets"] == False]
adata.obs.drop("predicted_doublets", axis=1, inplace=True)  # To avoid a potential TypeError

SHAPES.append(["DOUBLET", adata.shape[0], "-"])  

SHAPES


### 8. Plot after QC

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], save="_filterQC_counts_plot.pdf")  

sns.histplot(adata.obs["pct_counts_mt"])
plt.savefig("filterQC1_plot.pdf")

sns.displot(adata.obs["total_counts"], bins=100, kde=False)
plt.savefig("filterQC2_displot_total_counts.pdf")

sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=True,save="_filterQC3_plot.pdf")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=True,save="_rfilterQC4_plot.pdf")

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",save="_filterQC5_plot.pdf")

sns.histplot(adata.obs["pct_counts_ribo"])
plt.savefig("filterQC6_plot.pdf")


### 9. Output filter adata


# Output
SHAPES.append(["FINAL", *adata.shape])

with open(f"{OUTPUT_PATH}/QC.log", "w") as fo:
    for i, j, k in SHAPES:
        fo.write(f"{i}\t{j}\t{k}\n")

with open(f"{OUTPUT_PATH}/genes.txt", "w") as fo:
    fo.write("\n".join(adata.var_names))

adata.write_h5ad(f"{OUTPUT_PATH}/adata_filter.h5ad")




