# -------------------------------
# Load libraries
# -------------------------------
library(limma)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)
library(data.table)

# -------------------------------
# Read TPM matrix
# -------------------------------
counts <- read.csv("counts_matched_to_clinical_maf.csv", stringsAsFactors = FALSE)
counts$gene_name <- make.unique(as.character(counts$gene_name))  # Fix duplicates
rownames(counts) <- counts$gene_name
counts$gene_name <- NULL

# Remove leading X and replace dots with dashes (if any)
colnames(counts) <- gsub("^X", "", colnames(counts))
colnames(counts) <- gsub("\\.", "-", colnames(counts))
colnames(counts) <- trimws(colnames(counts))

# -------------------------------
# Read metadata
# -------------------------------
meta <- read.csv("meta_count_maf.csv", stringsAsFactors = FALSE, check.names = FALSE)
rownames(meta) <- meta$Sample
meta$Sample <- NULL
meta <- meta[trimws(rownames(meta)), , drop = FALSE]  # trim spaces

# -------------------------------
# Align samples
# -------------------------------
common_samples <- intersect(colnames(counts), rownames(meta))
if(length(common_samples) == 0){
  stop("No matching samples found between counts and metadata!")
}

# Subset both counts and meta to common samples
counts <- counts[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

# Check alignment
if(!all(colnames(counts) == rownames(meta))){
  stop("Samples are still not aligned! Please check your metadata and counts.")
}

# -------------------------------
# Create design matrix
# -------------------------------
meta$condition <- factor(meta$condition)
design <- model.matrix(~0 + condition, data = meta)
colnames(design) <- levels(meta$condition)

# -------------------------------
# voom transformation
# -------------------------------
dge <- DGEList(counts = counts)
dge$counts <- counts + 1  # TPM + 1 to avoid log(0)
v <- voom(dge, design = design, plot = TRUE)

# -------------------------------
# Linear modeling
# -------------------------------
fit <- lmFit(v, design)

# Automatically detect conditions for contrast
conds <- levels(meta$condition)
if(length(conds) != 2){
  stop("Currently only supports two-group comparison. Found conditions: ", paste(conds, collapse = ", "))
}
contrast.matrix <- makeContrasts(Contrast = paste0(conds[2], "-", conds[1]), levels = design)


fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# -------------------------------
# DEG results
# -------------------------------
res <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
res$gene <- rownames(res)
write.csv(res, "limma_all_results_post_clinical_maf_filter.csv", row.names = FALSE)


# -------------------------------
# Upregulated genes only (logFC > 1, adj.P.Val < 0.05)
# -------------------------------
up_res <- res[res$logFC > 1 & res$adj.P.Val < 0.05, ]
write.csv(up_res, "limma_upregulated_genes_after_clinical_filter.csv", row.names = FALSE)

cat("Upregulated genes saved:", nrow(up_res), "\n")


# -------------------------------
# MA plot
# -------------------------------
plotMA(fit2, coef = 1, main = paste0(conds[2], " vs ", conds[1]))

# -------------------------------
# PCA plot
# -------------------------------
pca <- prcomp(t(v$E))  # samples in rows
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)
pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     condition = meta$condition)

ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw()

# -------------------------------
# Heatmap (top 50 DEGs)
# -------------------------------
topgenes <- rownames(res[order(res$adj.P.Val), ])[1:50]
mat <- v$E[topgenes, ]

pheatmap(mat,
         scale = "row",
         annotation_col = meta,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# -------------------------------
# Volcano plot
# -------------------------------
# Volcano plot saved as PDF
# -------------------------------
pdf("Volcano_Tumor_vs_Control.pdf", width = 8, height = 6)  # open PDF device
suppressWarnings(
  EnhancedVolcano(res,
                  lab = res$gene,
                  x = 'logFC',
                  y = 'adj.P.Val',
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  title = "Tumor vs Control",
                  subtitle = "limma-voom differential expression")
)
dev.off()  # close PDF device


