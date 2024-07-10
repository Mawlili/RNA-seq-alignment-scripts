library(DESeq2)
library(MASS)
library(dplyr)

load("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/salmon_data_concise.RData")
#filter out NA 
pfas_w1_cord <- pfas_w1[!is.na(pfas_w1$cord_PFBA), ]
pfas_w1_mat <- pfas_w1[!is.na(pfas_w1$mat_PFBA), ]
#filter for common sample id
pfas_w1_cord_J <- paste0("J", pfas_w1_cord$condition.ID)
pfas_w1_mat_J <- paste0("J", pfas_w1_mat$condition.ID)
normalized_counts_filtered_cord <- normalized_counts[, colnames(normalized_counts) %in% pfas_w1_cord_J ]
normalized_counts_filtered_mat <- normalized_counts[, colnames(normalized_counts) %in% pfas_w1_mat_J ]
#construct deseqdataset and perfor deseq analysis
cord_pfas_diff <- DESeqDataSetFromMatrix(countData = normalized_counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + cord_PFBA)
mat_pfas_diff <- DESeqDataSetFromMatrix(countData = normalized_counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFBA)
dds_cord_pfas <- DESeq(cord_pfas_diff)
dds_mat_pfas <- DESeq(mat_pfas_diff)
res_cord_pfas <- results(cord_pfas_diff)
res_mat_pfas <- results(mat_pfas_diff)
save(res_cord_pfas, res_mat_pfas, file = "/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/pfas_results.RData")
