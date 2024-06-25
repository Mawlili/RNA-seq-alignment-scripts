library(DESeq2)
library(tximeta)

gse <- summarizeToGene(se)
expression_data <- assay(gse)
expression_data_log <- log2(expression_data + 1)
pca_result <- prcomp(t(expression_data_log))
png(file = "pca_plot_7.png", width = 800, height = 600, res = 120)
plot(pca_result$x[,1], pca_result$x[,2],
     xlab = "PC1", ylab = "PC2",
     main = "PCA of Log2-transformed Expression Data")
dev.off()
