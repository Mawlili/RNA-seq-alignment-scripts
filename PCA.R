#library
library(DESeq2)
library(tximeta)
library(RUVSeq)
library(MASS)
library(ggplot2)

#selecting house keeping gene
gse <- summarizeToGene(se)
counts <- round(assay(gse))
condition.GA <- coldata$condition.GA
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition.GA)
dds <- DESeq(dds)
res <- results(dds)
housekeeping_genes <- rownames(res[which(res$padj > 0.05), ])

#remove variance
set <- RUVg(counts, housekeeping_genes, k=1)
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
  xlim(-150, 150) +
  ylim(-50, 50)
  theme_minimal()
 ggsave("/rsrch5/home/canbio/whwu1/pca_plot_k1_2.png", plot = pca_plot, width = 10, height = 8, units = "in", dpi = 300) 

#DESeqDataSetFromMatrix w1+sex+...
diff_exp_transcript <- DESeqDataSetFromMatrix(countData = normalized_counts, colData = coldata, design = ~ W_1 + condition.sex + condition.ethnicity + condition.GROUP + condition.GA)
diff_exp_isoform <- DESeqDataSetFromMatrix(countData = se, colData = coldata, design = ~ W_1 + condition.sex + condition.ethnicity + condition.GROUP + condition.GA)


#gse and se

