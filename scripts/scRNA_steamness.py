import pandas as pd
import scanpy as sc
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import seaborn as sns
from cytotrace2_py.cytotrace2_py import cytotrace2
import os

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=600, facecolor="white")

# =====================================================
# LOAD DATA
# =====================================================
adata = sc.read("your_data.h5ad")

meta = pd.read_csv("meta.tsv", sep="\t", index_col=0)
common = adata.obs_names.intersection(meta.index)
adata = adata[common].copy()
adata.obs = meta.loc[common].copy()

# =====================================================
# QC FILTERING
# =====================================================
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True)

adata = adata[(adata.obs.pct_counts_mt < 25) & (adata.obs.pct_counts_ribo < 25), :]
sc.pp.filter_cells(adata, min_genes=250)
sc.pp.filter_genes(adata, min_cells=3)

# =====================================================
# NORMALIZATION
# =====================================================
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# =====================================================
# UMAP BEFORE STEMNESS
# =====================================================
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

sc.pl.umap(adata, color="Cell_subtype", frameon=False, save="_before_cytotrace2.png")

# =====================================================
# SUBSET: CANCER + NORMAL EPITHELIAL
# =====================================================
cancer_epi = adata.obs["Cell_subtype_original"].isin(["CMS1", "CMS2", "CMS3", "CMS4"])
normal_epi = (adata.obs["Cell_type"] == "Epithelial cells") & (adata.obs["Class"] == "Normal")

adata_epi = adata[cancer_epi | normal_epi].copy()

adata_epi.obs["Condition"] = "Normal"
adata_epi.obs.loc[cancer_epi, "Condition"] = "Cancer"
adata_epi.obs["Condition"] = adata_epi.obs["Condition"].astype("category")

print(adata_epi.obs["Condition"].value_counts())

# =====================================================
# WRITE REQUIRED FILES FOR CytoTRACE2
# =====================================================
print("\nWriting CytoTRACE2 input files...")

X = adata_epi.X
if sparse.issparse(X):
    X = X.toarray()

# ---- Expression file (genes × cells)
expr = pd.DataFrame(X.T, index=adata_epi.var_names, columns=adata_epi.obs_names)
expr_file = "cytotrace2_expression.tsv"
expr.to_csv(expr_file, sep="\t")

# ---- Annotation file (cell → group)
anno = adata_epi.obs[["Condition"]].copy()
anno.columns = ["Group"]
anno_file = "cytotrace2_annotation.tsv"
anno.to_csv(anno_file, sep="\t")

# =====================================================
# RUN CytoTRACE2
# =====================================================
print("\nRunning CytoTRACE2...")

cytotrace2(
    expr_file,
    annotation_path=anno_file,
    species="human"
)

# =====================================================
#Load CytoTRACE2 - inspect & correct columns
# =====================================================
print("\nLoading & inspecting CytoTRACE2 results...")

potency = pd.read_csv("cytotrace2_results.txt", sep="\t", index_col=0)

print("Potency DataFrame shape:", potency.shape)
print("\nFirst 3 columns:", potency.columns[:3].tolist())
print("\nColumn names (first 10):")
print(potency.columns[:10].tolist())
print("\nSample rows:")
print(potency.head(3))

# Correct column names (camelCase, no underscore based on snippet)
score_col = 'CytoTRACE2Score'  # Adjust if different, e.g., 'CytoTRACE2_Score'
potency_col = 'CytoTRACE2Potency'

if score_col not in potency.columns:
    print(f"\nAvailable columns containing 'Score':")
    score_cols = [col for col in potency.columns if 'Score' in col]
    print(score_cols)
    score_col = score_cols[0] if score_cols else potency.columns[1]  # Fallback to 2nd col
    print(f"Using '{score_col}' for score")

if potency_col not in potency.columns:
    print(f"\nAvailable columns containing 'Potency':")
    potency_cols = [col for col in potency.columns if 'Potency' in col]
    print(potency_cols)
    potency_col = potency_cols[0] if potency_cols else None
    print(f"Using '{potency_col}' for potency")

# Align & assign (handle NaNs from unmatched cells)
adata_epi.obs['CytoTRACE2_Score'] = potency.reindex(adata_epi.obs_names)[score_col]
adata_epi.obs['CytoTRACE2_Potency'] = potency.reindex(adata_epi.obs_names)[potency_col]

print("\nAssigned scores to", (~adata_epi.obs['CytoTRACE2_Score'].isna()).sum(), "cells")
print("Unique potencies:", adata_epi.obs['CytoTRACE2_Potency'].dropna().unique())

# Continue with stem definition & DE (same as before)
potency_levels = ['Unipotent', 'Oligopotent', 'Multipotent']  # Non-differentiated stem/progenitors[file:25]
cancer_mask = adata_epi.obs['Condition'] == 'Cancer'
stem_cells = adata_epi.obs['CytoTRACE2_Potency'].isin(potency_levels) & cancer_mask
normal_cells = adata_epi.obs['Condition'] == 'Normal'

print(f"\nCancer stem (non-Differentiated): {stem_cells.sum()}")
print(f"Normal epithelial: {normal_cells.sum()}")

adata_de = adata_epi[stem_cells | normal_cells].copy()
adata_de.obs['Group'] = 'Normal'
adata_de.obs.loc[stem_cells, 'Group'] = 'Cancer_Stem'
adata_de.obs['Group'] = adata_de.obs['Group'].astype('category')

print(adata_de.obs['Group'].value_counts())

if adata_de.shape[0] > 100:
    sc.pp.highly_variable_genes(adata_de, n_top_genes=2000)
    sc.pp.regress_out(adata_de, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata_de, max_value=10)
    sc.tl.pca(adata_de, n_comps=30)
    sc.pp.neighbors(adata_de, n_neighbors=10, n_pcs=30)
    sc.tl.umap(adata_de)

    sc.tl.rank_genes_groups(adata_de, 'Group', method='wilcoxon')

    result = sc.get.rank_genes_groups_df(adata_de, group='Cancer_Stem', pval_cutoff=0.05)
    result.to_csv('DE_CancerStem_vs_NormalEpi.csv')
    print("\nTop DE saved to DE_CancerStem_vs_NormalEpi.csv")

    sc.pl.rank_genes_groups(adata_de, n_genes=20, save='_CancerStem_DE.pdf')
    sc.pl.umap(adata_de, color=['Group', 'CytoTRACE2_Potency'], save='_CancerStem_UMAP.pdf')

print("\nFixed & done!")
")

# =====================================================
#  VOLCANO PLOT
# =====================================================
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv("DE_CancerStem_vs_NormalEpi.csv")

# Fix Scanpy column names if needed
if "logfoldchanges" in df.columns:
    df["log2FC"] = df["logfoldchanges"]
if "pvals_adj" in df.columns:
    df["padj"] = df["pvals_adj"]

fc_thresh = 1
p_thresh = 0.05

df["neglog10_padj"] = -np.log10(df["padj"])

# Categorize
df["Category"] = "NS"
df.loc[(df["log2FC"] > fc_thresh) & (df["padj"] < p_thresh), "Category"] = "Up"
df.loc[(df["log2FC"] < -fc_thresh) & (df["padj"] < p_thresh), "Category"] = "Down"

# ---------- PLOT ----------
plt.figure(figsize=(10,8))

# Background non-significant (soft light grey)
plt.scatter(
    df[df["Category"]=="NS"]["log2FC"],
    df[df["Category"]=="NS"]["neglog10_padj"],
    color="#DADADA",
    s=20,
    alpha=0.5
)

# Downregulated (soft lavender-blue)
plt.scatter(
    df[df["Category"]=="Down"]["log2FC"],
    df[df["Category"]=="Down"]["neglog10_padj"],
    color="#6A9AE2",
    s=55,
    alpha=0.85
)

# Upregulated (soft rose-magenta)
plt.scatter(
    df[df["Category"]=="Up"]["log2FC"],
    df[df["Category"]=="Up"]["neglog10_padj"],
    color="#E76F92",
    s=55,
    alpha=0.85
)

# Threshold lines
plt.axvline(fc_thresh, linestyle="--", color="black", lw=1)
plt.axvline(-fc_thresh, linestyle="--", color="black", lw=1)
plt.axhline(-np.log10(p_thresh), linestyle="--", color="black", lw=1)

# Symmetric X axis
max_fc = np.max(np.abs(df["log2FC"]))
plt.xlim(-max_fc*1.1, max_fc*1.1)

# Clean theme
plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)

plt.xlabel("log2 Fold Change", fontsize=14)
plt.ylabel("-log10 Adjusted P-value", fontsize=14)
plt.title("Cancer Stem vs Normal Epithelial Cells", fontsize=16, weight="bold")

# Label top genes
top_genes = df[(df["Category"]!="NS")].nsmallest(10, "padj")

for _, row in top_genes.iterrows():
    plt.text(
        row["log2FC"],
        row["neglog10_padj"]+0.3,
        row["names"],
        fontsize=9,
        ha='center'
    )

plt.tight_layout()
plt.savefig("Volcano_CancerStem_vs_Normal_PREMIUM.png", dpi=600)
plt.show()

print("🔥 Premium volcano plot saved as Volcano_CancerStem_vs_Normal_PREMIUM.png")


# =====================================================
# 9. VISUALIZATION PANELS
# =====================================================
print("\nCreating premium visualizations...")

sc.settings.set_figure_params(dpi=600, facecolor="white", fontsize=12)

# -----------------------------------------------------
# Epithelial vs Other
# -----------------------------------------------------
adata.obs["Epithelial_status"] = "Other"
adata.obs.loc[adata.obs["Cell_type"] == "Epithelial cells", "Epithelial_status"] = "Epithelial"
adata.obs["Epithelial_status"] = adata.obs["Epithelial_status"].astype("category")

sc.pl.umap(
    adata,
    color="Epithelial_status",
    palette={"Epithelial": "#E64B35", "Other": "#EFEFEF"},
    size=0.7,
    add_outline=True,
    outline_width=(0.12, 0.03),
    frameon=False,
    title="Epithelial Cells",
    save="_01_epithelial_vs_other.png"
)

# -----------------------------------------------------
# CMS-like groups from Cell_subtype (if available)
# -----------------------------------------------------
adata.obs["CMS_group"] = adata.obs["Cell_subtype"].astype(str)
adata.obs.loc[~adata.obs["CMS_group"].isin(["CMS1", "CMS2", "CMS3", "CMS4"]), "CMS_group"] = "Other"
adata.obs["CMS_group"] = adata.obs["CMS_group"].astype("category")

cms_colors = {
    "CMS1": "#3B4CC0",
    "CMS2": "#E64B35",
    "CMS3": "#4DBBD5",
    "CMS4": "#00A087",
    "Other": "#F0F0F0"
}

sc.pl.umap(
    adata,
    color="CMS_group",
    palette=cms_colors,
    size=0.7,
    add_outline=True,
    outline_width=(0.12, 0.03),
    frameon=False,
    title="CMS-like Subtypes in Epithelial Cells",
    save="_02_cms_subtypes.png"
)

# -----------------------------------------------------
# Tumor vs Normal vs Other
# -----------------------------------------------------
adata.obs["Tumor_status"] = "Other"
adata.obs.loc[adata.obs["Class"] == "Tumor", "Tumor_status"] = "Tumor"
adata.obs.loc[adata.obs["Class"] == "Normal", "Tumor_status"] = "Other"
adata.obs["Tumor_status"] = adata.obs["Tumor_status"].astype("category")

sc.pl.umap(
    adata,
    color="Tumor_status",
    palette={"Tumor": "#DC0000", "Normal": "#0072B2", "Other": "#EFEFEF"},
    size=0.7,
    add_outline=True,
    frameon=False,
    title="Tumor vs Other Cells",
    save="_03_tumor_vs_other.png"
)

# -----------------------------------------------------
# All Cell Subtypes
# -----------------------------------------------------
sc.pl.umap(
    adata,
    color="Cell_subtype",
    palette="tab20",
    size=0.7,
    add_outline=True,
    outline_width=(0.12, 0.03),
    frameon=False,
    legend_loc="right margin",
    legend_fontsize=8,
    title="Global Cell Subtype Atlas",
    save="_04_all_subtypes.png"
)


# =====================================================
# PATIENT COMPOSITION ANALYSIS 
# =====================================================
print("\nCalculating patient distribution across cell subtypes...")

patient_col = "Patient"

composition = pd.crosstab(adata.obs["Cell_subtype"], adata.obs[patient_col])
composition_perc = composition.div(composition.sum(axis=1), axis=0) * 100

os.makedirs("figures", exist_ok=True)

# --- A. Stacked Bar Plot ---
plt.figure(figsize=(13,6))
composition_perc.plot(
    kind="bar",
    stacked=True,
    colormap="tab20",
    width=0.85,
    edgecolor="white",
    linewidth=0.4
)

plt.ylabel("Percentage of Cells")
plt.xlabel("Cell Subtype")
plt.title("Patient Composition Across Cell Subtypes")
plt.legend(title="Patient", bbox_to_anchor=(1.02,1), loc="upper left", frameon=False)
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("figures/patient_composition_bar_PREMIUM.png", dpi=600)
plt.show()

# --- B. Heatmap ---
plt.figure(figsize=(11,8))
sns.heatmap(
    composition,
    cmap="Blues",
    linewidths=0.3,
    linecolor="white",
    cbar_kws={"label": "Cell Count"}
)
plt.title("Cell Counts per Patient and Subtype")
plt.ylabel("Cell Subtype")
plt.xlabel("Patient")
plt.tight_layout()
plt.savefig("figures/patient_counts_heatmap_PREMIUM.png", dpi=600)
plt.show()

# --- C. UMAP by Patient ---
sc.pl.umap(
    adata,
    color=patient_col,
    palette="Set2",
    frameon=False,
    size=0.7,
    title="Patient Distribution Across UMAP",
    save="_05_umap_by_patient.png"
)
