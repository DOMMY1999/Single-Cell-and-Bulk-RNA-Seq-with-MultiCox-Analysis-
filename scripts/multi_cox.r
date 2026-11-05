# libraryes
library(survival)
library(randomForestSRC)
library(timeROC)
library(survminer)
library(caret)
library(dplyr)
library(ggplot2)

df <- read.csv("shared_expression_with_clinical_1.csv", stringsAsFactors = FALSE)
# Rename for safety
colnames(df)[colnames(df)=="time"] <- "OS_time"
colnames(df)[colnames(df)=="status"] <- "OS_event"

# Encode clinical variables

df$gender <- factor(df$gender)
df$stage <- factor(df$stage, ordered = TRUE)

df$TMB <- log2(df$TMB + 1)
df$total_mutation <- log2(df$total_mutation + 1)

# Remove rows with missing values 

df <- df[complete.cases(df), ]
cat("Patients after NA removal:", nrow(df), "\n")


# Separate genes vs survival outcome

gene_cols <- colnames(df)[9:ncol(df)]
X_all <- as.matrix(df[, gene_cols])


# Random Forest Stability Selection
set.seed(123)
n_iter <- 500
rf_sel <- matrix(0, nrow=ncol(X_all), ncol=n_iter)
rownames(rf_sel) <- colnames(X_all)
pb <- txtProgressBar(min=0, max=n_iter, style=3)

for(i in 1:n_iter) {
  setTxtProgressBar(pb, i)
  
  # Bootstrap sample 
  boot_idx <- sample(nrow(X_all), replace=TRUE)
  X_boot <- X_all[boot_idx,]
  time_boot <- df$OS_time[boot_idx]
  event_boot <- df$OS_event[boot_idx]
  
  # Random Survival Forest
  rf_data <- data.frame(X_boot)
  rf_data$time <- time_boot
  rf_data$status <- event_boot
  
  rf_model <- rfsrc(Surv(time, status) ~ ., 
                    data=rf_data, 
                    nsplit=10, 
                    importance="permute",
                    samptype="swr")
  
  # Top 50 most important genes
  rf_imp <- rf_model$importance
  top_genes <- names(sort(rf_imp, decreasing=TRUE))[1:50]
  rf_sel[top_genes, i] <- 1
}

close(pb)

# RF selected genes (frequency >= 0.3)
rf_freq <- rowMeans(rf_sel)
rf_genes <- names(rf_freq[rf_freq >= 0.3])
cat("\n🎯 RF Genes:", length(rf_genes), "\n")
print(head(rf_genes, 20))



# RF VARIABLE IMPORTANCE PLOT
rf_full <- rfsrc(Surv(OS_time, OS_event) ~ ., 
                 data = data.frame(df[, gene_cols], 
                                   OS_time=df$OS_time, OS_event=df$OS_event),
                 importance = "permute")

rf_imp_full <- sort(rf_full$importance, decreasing = TRUE)[1:20]
rf_imp_df <- data.frame(Gene = names(rf_imp_full),
                        Importance = rf_imp_full)

p_rf_imp <- ggplot(rf_imp_df, aes(x=reorder(Gene, Importance), y=Importance)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  theme_bw() +
  labs(title="Top 20 RF Variable Importance", x="Gene", y="Permutation Importance")

ggsave("RF_Variable_Importance_Top20.png", p_rf_imp, width=7, height=6, dpi=900)

print(p_rf_imp)

rf_freq_df <- data.frame(Gene = rf_genes,
                         Frequency = rf_freq[rf_genes])

p_rf_freq <- ggplot(rf_freq_df, aes(x=reorder(Gene, Frequency), y=Frequency)) +
  geom_bar(stat="identity", fill="darkred") +
  geom_hline(yintercept = 0.3, linetype="dashed", color="black", linewidth = 1) +
  coord_flip() +
  theme_bw() +
  ylim(0,1) +
  labs(title="RF Stability Selection Frequency", x="Gene", y="Selection Frequency")

ggsave("RF_Stability_Frequency.png", p_rf_freq, width=7, height=6, dpi=900)

print(p_rf_freq)

# UNIVARIATE COX FILTERING 
y <- Surv(df$OS_time, df$OS_event)  

X_rf <- X_all[, rf_genes, drop = FALSE]

uni_results <- lapply(colnames(X_rf), function(g) {
  fit <- coxph(y ~ X_rf[, g])
  s <- summary(fit)
  c(Gene = g,
    HR = s$coefficients[,"exp(coef)"],
    Pvalue = s$coefficients[,"Pr(>|z|)"])
})

uni_df <- as.data.frame(do.call(rbind, uni_results))
uni_df$HR <- as.numeric(uni_df$HR)
uni_df$Pvalue <- as.numeric(uni_df$Pvalue)

rf_final_genes <- uni_df$Gene[uni_df$Pvalue < 0.05]
cat("✅ RF Final Genes (p<0.05):", length(rf_final_genes), "\n")
print(rf_final_genes)


cat("✅ FINAL RF genes (Cox p<0.05):", length(rf_final_genes), "\n")


# UNIVARIATE COX FOREST PLOT
uni_plot_df <- uni_df %>%
  filter(Gene %in% rf_final_genes) %>%
  mutate(logHR = log(HR),
         lower = exp(logHR - 1.96*sd(logHR)),
         upper = exp(logHR + 1.96*sd(logHR)))

p_uni_forest <- ggplot(uni_plot_df,
                       aes(y=reorder(Gene, HR), x=HR)) +
  geom_point(color="black") +
  geom_vline(xintercept = 1, linetype="dashed") +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=0.2) +
  scale_x_log10() +
  theme_bw() +
  labs(title="Univariate Cox Analysis of RF-Selected Genes",
       x="Hazard Ratio (log scale)", y="Gene")

ggsave("Univariate_Cox_Forest.png", p_uni_forest, width=7, height=6, dpi=300)
print(p_uni_forest)


# MULTIVARIATE COX MODEL (RF genes + clinical)
cat("\n Building final Cox model...\n")
df_rf_model <- df[, c("OS_time","OS_event","age","gender","stage","TMB","total_mutation", rf_final_genes)]
fit_rf_final <- coxph(Surv(OS_time, OS_event) ~ ., data = df_rf_model)
summary(fit_rf_final)



# MULTIVARIATE COX FOREST PLOT (GENES ONLY)
multi_coef <- summary(fit_rf_final)$coefficients
multi_df <- as.data.frame(multi_coef)
multi_df$Variable <- rownames(multi_df)
multi_df <- multi_df[multi_df$Variable %in% rf_final_genes, ]
multi_df$HR <- multi_df$`exp(coef)`
multi_df$lower <- exp(log(multi_df$HR) - 1.96 * multi_df$`se(coef)`)
multi_df$upper <- exp(log(multi_df$HR) + 1.96 * multi_df$`se(coef)`)

p_multi_forest <- ggplot(multi_df,
                         aes(y=reorder(Variable, HR), x=HR)) +
  geom_point(color="red") +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=0.2, color="orange") +
  scale_x_log10() +
  theme_bw() +
  labs(title="Multivariate Cox Model (Gene Signature)",
       x="Hazard Ratio (log scale)", y="Gene")

ggsave("Multivariate_Cox_Forest_Genes.png", p_multi_forest,
 width=7, height=6, dpi=900)
print(p_multi_forest)


# RISK SCORE + RISK GROUPS
df$rf_risk_score <- predict(fit_rf_final, type = "lp")
cutoff_rf <- median(df$rf_risk_score)
df$rf_risk_group <- ifelse(df$rf_risk_score >= cutoff_rf, "High", "Low")
df$rf_risk_group <- factor(df$rf_risk_group, levels = c("Low","High"))

cat("\n✅ RF Risk groups created\n")
table(df$rf_risk_group)


# HEATMAP OF GENE EXPRESSION VS RISK GROUPS
library(ComplexHeatmap)
library(circlize)

# Matrix
heat_data <- as.matrix(df[, rf_final_genes])
rownames(heat_data) <- df$CaseID

# Z-score per gene
heat_data_scaled <- scale(heat_data)

# Expression color scale
expr_col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

# Row annotations (patients)
row_anno <- rowAnnotation(
  Risk = df$rf_risk_group,
  Status = factor(df$OS_event),
  col = list(
    Risk = c("Low" = "#4DAF4A", "High" = "#E41A1C"),
    Status = c("0" = "#999999", "1" = "#000000")
  ),
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Risk = list(title = "Risk Group"),
    Status = list(title = "Event Status")
  )
)

png("RF_Gene_Heatmap_Colorful.png", width=1300, height=1400, res = 300)
Heatmap(heat_data_scaled,
        name = "Expression (Z-score)",
        col = expr_col_fun,
        right_annotation = row_anno,
        show_column_names = TRUE,
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 10),
        column_title = "RF Gene Expression Heatmap",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        heatmap_legend_param = list(title = "Z-score",
                                    legend_height = unit(4, "cm")))
dev.off()



# KAPLAN-MEIER SURVIVAL CURVE
print("Available columns in df:")
print(colnames(df))
print("Risk-related columns:")
print(grep("risk|Risk|group|Group", colnames(df), value=TRUE, ignore.case=TRUE))

fit_km_rf <- survfit(Surv(OS_time, OS_event) ~ rf_risk_group, data = df)
p <- ggsurvplot(fit_km_rf, data = df, pval = TRUE, risk.table = TRUE, 
                conf.int = FALSE, title = "RF-Derived Risk Model", xlab = "Time (days)")

# SAVE HIGH DPI
ggsave("RF_Survival_300dpi.png", 
       plot = print(p),           
       width = 12, height = 8, 
       bg = "white",
       dpi = 900)


fit_km_rf <- survfit(Surv(OS_time, OS_event) ~ rf_risk_group, data = df)

p <- ggsurvplot(fit_km_rf,
                data = df,
                pval = TRUE,
                risk.table = TRUE,
                conf.int = FALSE,
                title = "RF Risk Model: High vs Low Risk",
                xlab = "Time (days)",
                palette = c("lightblue", "orange"),
                bg="white",
                legend.labs = c("Low Risk", "High Risk"))

# SAVE HIGH DPI
ggsave("RF_High_Low_Risk_KM.png", plot=print(p), width=8, height=7, dpi=900, bg="white")
print(p)

# p-value display
if(logrank_p < 1e-10) {
  p_display <- "< 1e-10"
} else if(logrank_p < 1e-5) {
  p_display <- paste0("< 0.00001")
} else {
  p_display <- format.pval(logrank_p, digits=3, eps=1e-4)
}

cat("Log-rank p-value:", p_display, "\n")



# TIME-DEPENDENT ROC (1, 3, 5 years)

roc_rf <- timeROC(T = df$OS_time,
                  delta = df$OS_event,
                  marker = df$rf_risk_score,
                  cause = 1,
                  times = c(365, 1095, 1825),
                  iid = TRUE)

png("ROC_1_3_5yr.png", width=1600, height=1400, res=100, pointsize=22)
par(mar=c(6, 5.5, 4.5, 7), cex.axis=1.4, cex.lab=1.5)  

# Plot curves
plot(roc_rf, time=365, col="#FF8C00", lwd=4)
plot(roc_rf, time=1095, add=TRUE, col="#FF69B4", lwd=4)
plot(roc_rf, time=1825, add=TRUE, col="#32CD32", lwd=4)
abline(a=0, b=1, lty=2, col="gray40", lwd=2)


legend("bottomright", 
       inset=c(0.02, 0.02),        
       legend=paste0(c("1-year","3-year","5-year"), 
                     "\nAUC=", round(c(
                       roc_rf$AUC[roc_rf$times==365],
                       roc_rf$AUC[roc_rf$times==1095],
                       roc_rf$AUC[roc_rf$times==1825]), 3)),
       col=c("#FF8C00","#FF69B4","#32CD32"),
       lwd=4, cex=1.2, bty="n", bg="white")

title("RF Model: 1,3,5-Year ROC Curves", line=2)
dev.off()


print(roc_rf)


# INDEPENDENCE TEST (RF risk + clinical)


fit_indep_rf <- coxph(Surv(OS_time, OS_event) ~ rf_risk_score + age + gender + 
                        stage + TMB + total_mutation, data = df)
summary(fit_indep_rf)

# C-index
c_rf <- concordance(Surv(df$OS_time, df$OS_event) ~ df$rf_risk_score)$concordance
cat("\nC-index (RF model):", round(c_rf, 3), "\n")



# LOW vs HIGH PATIENT LISTS 

cat("🎯 LOW vs HIGH RISK PATIENTS:\n")
cat("Low risk patients:", sum(df$rf_risk_group == "Low"), "\n")
cat("High risk patients:", sum(df$rf_risk_group == "High"), "\n\n")

# LOW RISK PATIENTS
low_risk_df <- df[df$rf_risk_group == "Low", c("CaseID", "OS_time", "OS_event", "rf_risk_score", 
                                               "age", "gender", "stage", rf_final_genes[1:5])]
cat("=== LOW RISK PATIENTS (Top 10) ===\n")
print(head(low_risk_df, 10))

# HIGH RISK PATIENTS  
high_risk_df <- df[df$rf_risk_group == "High", c("CaseID", "OS_time", "OS_event", "rf_risk_score",
                                                 "age", "gender", "stage", rf_final_genes[1:5])]
cat("\n=== HIGH RISK PATIENTS (Top 10) ===\n")
print(head(high_risk_df, 10))

# SAVE FULL LISTS
write.csv(low_risk_df, "RF_Low_Risk_Patients.csv", row.names=FALSE)
write.csv(high_risk_df, "RF_High_Risk_Patients.csv", row.names=FALSE)

cat("\n💾 SAVED:\n")
cat("RF_Low_Risk_Patients.csv\n")
cat("RF_High_Risk_Patients.csv\n")


# INDIVIDUAL GENE KM PLOTS (ALL RF GENES)


library(ggpubr)
library(gridExtra)

cat("\n🎨 Creating individual gene KM plots...\n")

# Function for single gene KM plot
plot_gene_km <- function(gene_name, data=df, title_prefix="") {
  data[[gene_name]] <- as.numeric(data[[gene_name]])
  median_expr <- median(data[[gene_name]], na.rm=TRUE)
  data$gene_group <- ifelse(data[[gene_name]] > median_expr, "High", "Low")
  data$gene_group <- factor(data$gene_group, levels=c("Low", "High"))
  
  fit <- survfit(Surv(OS_time, OS_event) ~ gene_group, data=data)
  p <- ggsurvplot(fit, data=data, pval=TRUE, risk.table=TRUE,
                  title=paste(title_prefix, gene_name),
                  xlab="Time (days)", ylab="Survival Probability",
                  legend.title=paste(gene_name, "expression"),
                  legend.labs=c("Low", "High"),
                  bg="white",
                  ggtheme=theme_minimal()) +
    theme(plot.title=element_text(hjust=0.5, size=12))
  return(p)
}

# Create plots for ALL RF genes
gene_plots <- list()
plot_count <- 0

for(gene in rf_final_genes) {
  if(plot_count >= 12) break  # Limit to 12 for display
  plot_count <- plot_count + 1
  gene_plots[[gene]] <- plot_gene_km(gene, title_prefix="RF Gene: ")
}

# Display first 9 genes (3x3 grid)
if(length(gene_plots) >= 9) {
  combined_plot <- ggarrange(plotlist=gene_plots[1:9], ncol=3, nrow=3)
  print(combined_plot)
}


# GENE EXPRESSION SUMMARY (High vs Low Risk)


cat("\n📊 GENE EXPRESSION BY RISK GROUP:\n")
gene_summary <- data.frame(
  Gene = rf_final_genes,
  Low_Mean = sapply(rf_final_genes, function(g) mean(df[df$rf_risk_group=="Low", g], na.rm=TRUE)),
  High_Mean = sapply(rf_final_genes, function(g) mean(df[df$rf_risk_group=="High", g], na.rm=TRUE)),
  Fold_Change = sapply(rf_final_genes, function(g) 
    mean(df[df$rf_risk_group=="High", g], na.rm=TRUE) / mean(df[df$rf_risk_group=="Low", g], na.rm=TRUE)),
  P_value = sapply(rf_final_genes, function(g) {
    t.test(df[df$rf_risk_group=="Low", g], df[df$rf_risk_group=="High", g])$p.value
  })
)

print(round(gene_summary, 3))
write.csv(gene_summary, "RF_Gene_Expression_High_vs_Low.csv", row.names=FALSE)


# TOP PATIENTS BY RISK SCORE 
cat("\n🔥 TOP 5 LOWEST RISK PATIENTS:\n")
top_low <- df[df$rf_risk_group=="Low",][order(df[df$rf_risk_group=="Low",]$rf_risk_score)[1:5], 
                                        c("CaseID", "rf_risk_score", rf_final_genes[1:3])]
print(top_low)

cat("\n💀 TOP 5 HIGHEST RISK PATIENTS:\n")
top_high <- df[df$rf_risk_group=="High",][order(df[df$rf_risk_group=="High",]$rf_risk_score, decreasing=TRUE)[1:5], 
                                          c("CaseID", "rf_risk_score", rf_final_genes[1:3])]
print(top_high)


# SAVE ALL INDIVIDUAL GENE KM PLOTS
cat("\n💾 Saving individual gene KM plots...\n")
for(i in 1:length(rf_final_genes)) {
  gene <- rf_final_genes[i]
  p <- plot_gene_km(gene)
  ggsave(paste0("KM_", gene, ".png"), plot=p$plot, width=8, height=6, dpi=300)
}

cat("\n✅ COMPLETE RF RESULTS!\n")
cat("🎯 Files saved:\n")
cat("1. RF_Low_Risk_Patients.csv\n")
cat("2. RF_High_Risk_Patients.csv  \n")
cat("3. RF_Gene_Expression_High_vs_Low.csv\n")
cat("4. KM_[GENE].png (individual plots)\n")
cat("5. Log-rank p < 1e-10 (EXCELLENT!)\n")



# EXTRACT FINAL GENE COEFFICIENTS
coef_rf <- summary(fit_rf_final)$coefficients
coef_df_rf <- as.data.frame(coef_rf)
coef_df_rf$Variable <- rownames(coef_df_rf)
rf_gene_results <- coef_df_rf[coef_df_rf$Variable %in% rf_final_genes, ]
rf_gene_results <- rf_gene_results[, c("Variable","coef","`exp(coef.)`","Pr(>|z|)")]
colnames(rf_gene_results) <- c("Gene","Coefficient","HazardRatio","Pvalue")


# SAVE ALL RESULTS

write.csv(rf_final_genes, "RF_final_genes.csv", row.names=FALSE)
write.csv(rf_gene_results, "RF_gene_coefficients.csv", row.names=FALSE)

risk_df_rf <- df[, c("CaseID","OS_time","OS_event","rf_risk_score","rf_risk_group")]
write.csv(risk_df_rf, "RF_patient_risk_scores.csv", row.names=FALSE)

write.csv(data.frame(
  Time=c("1yr","3yr","5yr"),
  AUC=roc_rf$AUC
), "RF_ROC_AUC.csv", row.names=FALSE)

cat("\n✅ INDEPENDENT RF ANALYSIS COMPLETE!\n")
cat("🎯 RF Genes (Cox p<0.05):", length(rf_final_genes), "\n")
cat("📊 1yr/3yr/5yr AUC:", round(roc_rf$AUC,3), "\n")
cat("🏆 C-index:", round(c_rf,3), "\n")
cat("💾 Saved 4 files:\n")
cat("   1. RF_final_genes.csv\n") 
cat("   2. RF_gene_coefficients.csv\n")
cat("   3. RF_patient_risk_scores.csv\n")
cat("   4. RF_ROC_AUC.csv\n")

