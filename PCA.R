#library
library(DESeq2)
library(tximeta)
library(RUVSeq)
library(MASS)
library(ggplot2)
library(EnhancedVolcano)
library(readxl)
library(dplyr)
library(writexl)

#selecting house keeping gene
gse <- summarizeToGene(se)
counts_se <- round(assay(se))
counts <- round(assay(gse))
condition.GA <- coldata$condition.GA
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition.GA)
dds <- DESeq(dds)
res <- counts(dds, normalized = TRUE)
gene_variance <- apply(res, 1, var)
low_variance_genes <- names(sort(gene_variance)[1:100])

#housekeeping_genes <- rownames(res[which(res$padj > 0.05), ])

#remove variance
set <- RUVg(counts, low_variance_genes, k=1)
normalized_counts <- set$normalizedCounts

#vst normalization
dds_normalized <- DESeqDataSetFromMatrix(countData = normalized_counts, colData = coldata, design = ~ condition.GA)
vst_data <- vst(dds_normalized)

#log normalization and pca 
pca_data <- prcomp(t(assay(vst_data)))
pca_components <- data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2])
pca_plot_data <- cbind(pca_components, coldata)
pca_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color =  condition.sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(100 * summary(pca)$importance[2, 1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * summary(pca)$importance[2, 2], 1), "% variance")) +
  ggtitle("PCA Plot") +
 xlim(-20, 25) +
  ylim(-100, 100)
  theme_minimal()
 ggsave("pca_var.png", plot = pca_plot, width = 10, height = 8, units = "in", dpi = 300) 

#DESeqDataSetFromMatrix w1+sex+...
W_1 <- set$W
coldata_w1 <- cbind(coldata, D = W_1)
diff_exp_transcript <- DESeqDataSetFromMatrix(countData = normalized_counts, colData = coldata_w1, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA)
diff_exp_isoform <- DESeqDataSetFromMatrix(countData = counts_se, colData = coldata_w1, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA)
dds_transcript <- DESeq(diff_exp_transcript)
dds_isoform <- DESeq(diff_exp_isoform)
res_transcript <- results(dds_transcript)
res_isoform <- results(dds_isoform)


#pfas analysis
pfas <- read_excel("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/pfas.xlsx")
pfas <- pfas %>%
  rename(condition.SubjectID = SubjectID)
pfas_w1 <- inner_join(coldata_w1, pfas, by = "condition.SubjectID")
#pfba gene family for transcript

normalized_counts_f <- dds[rowSums(counts(dds)) > 1, ] 
cord_pfas_diff <- DESeqDataSetFromMatrix(countData = counts, colData = pfas_w1, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex)
mat_pfas_diff <- DESeqDataSetFromMatrix(countData = counts, colData = pfas_w1, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFBA)
dds_cord_pfas <- DESeq(cord_pfas_diff)
dds_mat_pfas <- DESeq(mat_pfas_diff)
res_cord_pfas <- results(cord_pfas_diff)
res_mat_pfas <- results(mat_pfas_diff)
