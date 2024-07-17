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
normalized_counts_filtered_cord <- normalized_counts[, colnames(counts) %in% pfas_w1_cord_J ]
normalized_counts_filtered_mat <- normalized_counts[, colnames(normalized_counts) %in% pfas_w1_mat_J ]

counts_filtered_cord <- counts[, colnames(counts) %in% pfas_w1_cord_J ]
counts_filtered_mat <- counts[, colnames(counts) %in% pfas_w1_mat_J ]
#construct deseqdataset and perfor deseq analysis

#cord PFBA
cord_pfas_diff1 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFBA)
cord_pfas_diff1 <- estimateSizeFactors(cord_pfas_diff1)
idx1 <- rowSums(counts(cord_pfas_diff1, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff1 <- cord_pfas_diff1[idx1,]
dds_cord_pfas1 <- DESeq(cord_pfas_diff1)
res_cord_pfas1 <- results(dds_cord_pfas1)

#cord PFOA
cord_pfas_diff2 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFOA)
cord_pfas_diff2 <- estimateSizeFactors(cord_pfas_diff2)
idx2 <- rowSums(counts(cord_pfas_diff2, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff2 <- cord_pfas_diff2[idx2,]
dds_cord_pfas2 <- DESeq(cord_pfas_diff2)
res_cord_pfas2 <- results(dds_cord_pfas2)

#cord PFNA
cord_pfas_diff3 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFNA)
cord_pfas_diff3 <- estimateSizeFactors(cord_pfas_diff3)
idx3 <- rowSums(counts(cord_pfas_diff3, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff3 <- cord_pfas_diff3[idx3,]
dds_cord_pfas3 <- DESeq(cord_pfas_diff3)
res_cord_pfas3 <- results(dds_cord_pfas3)

#cord PFDA
cord_pfas_diff4 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFDA)
cord_pfas_diff4 <- estimateSizeFactors(cord_pfas_diff4)
idx4 <- rowSums(counts(cord_pfas_diff4, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff4 <- cord_pfas_diff4[idx4,]
dds_cord_pfas4 <- DESeq(cord_pfas_diff4)
res_cord_pfas4 <- results(dds_cord_pfas4)

#cord PFUnDA
cord_pfas_diff5 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFUnDA)
cord_pfas_diff5 <- estimateSizeFactors(cord_pfas_diff5)
idx5 <- rowSums(counts(cord_pfas_diff5, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff5 <- cord_pfas_diff5[idx5,]
dds_cord_pfas5 <- DESeq(cord_pfas_diff5)
res_cord_pfas5 <- results(dds_cord_pfas5)

#cord PFBS
cord_pfas_diff6 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFBS)
cord_pfas_diff6 <- estimateSizeFactors(cord_pfas_diff6)
idx6 <- rowSums(counts(cord_pfas_diff6, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff6 <- cord_pfas_diff6[idx6,]
dds_cord_pfas6 <- DESeq(cord_pfas_diff6)
res_cord_pfas6 <- results(dds_cord_pfas6)

#cord PFHxS
cord_pfas_diff7 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFHxS)
cord_pfas_diff7 <- estimateSizeFactors(cord_pfas_diff7)
idx7 <- rowSums(counts(cord_pfas_diff7, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff7 <- cord_pfas_diff7[idx7,]
dds_cord_pfas7 <- DESeq(cord_pfas_diff7)
res_cord_pfas7 <- results(dds_cord_pfas7)

#mat PFBA
mat_pfas_diff1 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFBA)
mat_pfas_diff1 <- estimateSizeFactors(mat_pfas_diff1)
idx_m1 <- rowSums(counts(mat_pfas_diff1, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff1 <- mat_pfas_diff1[idx_m1,]
dds_mat_pfas1 <- DESeq(mat_pfas_diff1)
res_mat_pfas1 <- results(dds_mat_pfas1)

#mat PFOA
mat_pfas_diff2 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFOA)
mat_pfas_diff2 <- estimateSizeFactors(mat_pfas_diff2)
idx_m2 <- rowSums(counts(mat_pfas_diff2, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff2 <- mat_pfas_diff2[idx_m2,]
dds_mat_pfas2 <- DESeq(mat_pfas_diff2)
res_mat_pfas2 <- results(dds_mat_pfas2)

#mat PFNA
mat_pfas_diff3 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFNA)
mat_pfas_diff3 <- estimateSizeFactors(mat_pfas_diff3)
idx_m3 <- rowSums(counts(mat_pfas_diff3, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff3 <- mat_pfas_diff3[idx_m3,]
dds_mat_pfas3 <- DESeq(mat_pfas_diff3)
res_mat_pfas3 <- results(dds_mat_pfas3)

#mat PFDA
mat_pfas_diff4 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFDA)
mat_pfas_diff4 <- estimateSizeFactors(mat_pfas_diff4)
idx_m4 <- rowSums(counts(mat_pfas_diff4, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff4 <- mat_pfas_diff4[idx_m4,]
dds_mat_pfas4 <- DESeq(mat_pfas_diff4)
res_mat_pfas4 <- results(dds_mat_pfas4)

#mat PFUnDA
mat_pfas_diff5 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFUnDA)
mat_pfas_diff5 <- estimateSizeFactors(mat_pfas_diff5)
idx_m5 <- rowSums(counts(mat_pfas_diff5, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff5 <- mat_pfas_diff5[idx_m5,]
dds_mat_pfas5 <- DESeq(mat_pfas_diff5)
res_mat_pfas5 <- results(dds_mat_pfas5)

#mat PFBS
mat_pfas_diff6 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFBS)
mat_pfas_diff6 <- estimateSizeFactors(mat_pfas_diff6)
idx_m6 <- rowSums(counts(mat_pfas_diff6, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff6 <- mat_pfas_diff6[idx_m6,]
dds_mat_pfas6 <- DESeq(mat_pfas_diff6)
res_mat_pfas6 <- results(dds_mat_pfas6)

#mat PFHxS
mat_pfas_diff7 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFHxS)
mat_pfas_diff7 <- estimateSizeFactors(mat_pfas_diff7)
idx_m7 <- rowSums(counts(mat_pfas_diff7, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff7 <- mat_pfas_diff7[idx_m7,]
dds_mat_pfas7 <- DESeq(mat_pfas_diff7)
res_mat_pfas7 <- results(dds_mat_pfas7)

save(res_cord_pfas1, res_cord_pfas2, res_cord_pfas3, res_cord_pfas4, res_cord_pfas5, res_cord_pfas6, res_cord_pfas7, res_mat_pfas1, res_mat_pfas2, res_mat_pfas3, res_mat_pfas4, res_mat_pfas5, res_mat_pfas6, res_mat_pfas7, file = "/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/pfas_results_716.RData")
