library(DESeq2)
library(tximeta)

gse <- summarizeToGene(se)
expression_data <- assay(gse)
expression_data_log <- log2(expression_data + 1)
pca_plot <- ggplot(coldata, aes(x = PC1, y = PC2, color =  condition.mother_highest_education)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(100 * summary(pca)$importance[2, 1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * summary(pca)$importance[2, 2], 1), "% variance")) +
  ggtitle("PCA Plot") +
  xlim(-150, 150) +
  ylim(-200, 250)
  theme_minimal()
ggsave("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/pca_plot_10.png", plot = pca_plot, width = 10, height = 8, units = "in", dpi = 300)

