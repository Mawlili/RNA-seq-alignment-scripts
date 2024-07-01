#library
library(DESeq2)
library(tximeta)
library(RUVSeq)
library(MASS)

#selecting house keeping gene
gse <- summarizeToGene(se)
counts <- round(assay(gse))
condition.GA <- coldata$condition.GA
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition.GA)
dds <- DESeq(dds)
res <- results(dds)
housekeeping_genes <- rownames(res[which(res$padj > 0.05), ])
set <- RUVg(expression_data, housekeeping_genes, k=1)
#log normalization and pca 
pca_data <- prcomp(t(assay(rlog_data)))
pca_plot <- ggplot(coldata, aes(x = PC1, y = PC2, color =  condition.mother_highest_education)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(100 * summary(pca)$importance[2, 1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * summary(pca)$importance[2, 2], 1), "% variance")) +
  ggtitle("PCA Plot") +
  xlim(-150, 150) +
  ylim(-200, 250)
  theme_minimal()
ggsave("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/pca_plot_10.png", plot = pca_plot, width = 10, height = 8, units = "in", dpi = 300)

