Single-Cell and Bulk RNA-Seq MultiCox Analysis to Identify Cancer Stem Cell Specific Prognostic Genes associated with metastasis in CRC

Abstract: Metastatic progression and high recurrence rates remain major challenge of the colorectal cancer
(CRC), which is mainly driven my complex Tumor Microenvironment (TME) and high
heterogenic nature. Cancer stem cells (CSCs) are the subpopulation of TME characterized with
self-renewal, and pluripotent nature, they drive the tumor towards metastasis through cellular
plasticity and signal pathways interaction among TME. This study aims to screen CSC specific
biomarkers with metastatic risk, which can help guide patient risk stratification and improve
tumor prognosis through a multimodal approach that combines scRNA, bulk RNA, and
mutational profiling, CSC specific tumor cells within tumor clusters were identified using the
Cytotrace2 stemness scoring algorithm. Differentially expressed genes from scRNA and bulk
RNA were integrated with clinical features. These integrative clinical features were analyzed
using univariate Cox regression to identify gene expression signatures strongly associated with
patient status. Subsequently, multivariable Cox regression was performed to determine which
genes remain significant even after adjusting for clinical features. INA, UCN, FOXS1,
SEMA3E, UST, and EVX1 emerged as universal genes independently associated with poor
prognosis. Further in vitro validation confirmed that SEMA3E and FOXS1 were significantly
overexpressed in the cancer stem cells, and literature support suggests that SEMA3E and FOXS1
are strongly associated with metastasis in CRC.

Study objectives
- Characterize epithelial/tumor cell populations in CRC single-cell data.
- Identify genes differentially expressed between cancer and normal epithelial populations.
- Explore stemness-associated expression patterns.
- Integrate expression and clinical variables for survival analysis.
- Identify candidate prognostic genes and evaluate risk stratification using survival analyses, univarate and multivariate-cox models.
Repository contents
File	Description	Language
scRNA_steamness.py	scRNA-seq preprocessing, QC, normalization, PCA/UMAP, epithelial subset analysis, differential expression, stemness analysis, and visualizations	Python
multi_cox.r	Random Survival Forest selection, univariate and multivariable Cox analysis, risk scoring, Kaplan–Meier analysis, time-dependent ROC, and result export	R


Workflow
1. Single-cell analysis
The scRNA_stemcell.py script:
1. Reads an AnnData object and cell metadata.
2. Intersects cell identifiers and attaches metadata.
3. Groups CMS1–CMS4 labels under an epithelial category while retaining original subtype labels.
4. Calculates mitochondrial and ribosomal QC metrics and filters cells.
5. Filters low-quality cells and lowly expressed genes.
6. Normalizes counts and applies log transformation.
7. Performs PCA, neighborhood graph construction, and UMAP.
8. Subsets cancer and normal epithelial cells.
9. Performs cancer VS normal epithelial differential expression using the Wilcoxon method.
10. Generates differential expression plots and a heatmap.
11. Calculates a stmeness cluters using  CytoTRACE2 algorithm


Survival and prognostic gene analysis
The multi_cox.r script:
1. Loads expression and clinical survival data.
2. Encodes clinical variables and transforms mutation-related variables.
3. Removes rows containing missing values.
4. Performs bootstrap based Random Survival Forest feature-selection iterations.
5. Retains genes based on selection frequency.
6. Performs univariate Cox filtering.
7. Fits a multivariable Cox model incorporating selected genes and clinical variables.
8. Calculates a linear predictor risk score and divides patients into high- and low-risk groups.
9. Generates gene-importance plots, heatmaps, Kaplan–Meier curves, and time-dependent ROC curves.
10. Exports selected genes, coefficients, risk scores, and related analysis tables.

