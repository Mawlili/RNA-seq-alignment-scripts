library(DESeq2)
library(MASS)
library(dplyr)
library(readxl)

load("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/salmon_data_concise.RData")
#filter out NA 
pfas_w1 <- read_xlsx("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/pfas_w1_adjusted.xlsx")
pfas_w1_cord <- pfas_w1[!is.na(pfas_w1$cord_PFBA), ]
pfas_w1_mat <- pfas_w1[!is.na(pfas_w1$mat_PFBA), ]
pfas_w1_both <- pfas_w1_mat[!is.na(pfas_w1_mat$cord_PFBA), ]
#filter for common sample id
pfas_w1_cord_J <- paste0("J", pfas_w1_cord$condition.ID)
pfas_w1_mat_J <- paste0("J", pfas_w1_mat$condition.ID)
pfas_w1_both_J <- paste0("J", pfas_w1_both$condition.ID)

counts_filtered_cord <- counts[, colnames(counts) %in% pfas_w1_cord_J ]
counts_filtered_mat <- counts[, colnames(counts) %in% pfas_w1_mat_J ]
counts_filtered_both <- counts[, colnames(counts) %in% pfas_w1_both_J ]

#construct deseqdataset and perform deseq analysis

#PFBA
c_pfas_diff1 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + mat_PFBA + cord_PFBA)
c_pfas_diff1 <- estimateSizeFactors(c_pfas_diff1)
idxc1 <- rowSums(counts(c_pfas_diff1, normalized=TRUE) >= 5 ) >= 3
c_pfas_diff1 <- c_pfas_diff1[idxc1,]
dds_c_pfas1 <- DESeq(c_pfas_diff1)
res_c_pfas1 <- results(dds_c_pfas1)

#PFOA
c_pfas_diff2 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + mat_PFOA + cord_PFOA)
c_pfas_diff2 <- estimateSizeFactors(c_pfas_diff2)
idxc2 <- rowSums(counts(c_pfas_diff2, normalized=TRUE) >= 5 ) >= 3
c_pfas_diff2 <- c_pfas_diff2[idxc2,]
dds_c_pfas2 <- DESeq(c_pfas_diff2)
res_c_pfas2 <- results(dds_c_pfas2)

#PFNA
c_pfas_diff3 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + mat_PFNA + cord_PFNA)
c_pfas_diff3 <- estimateSizeFactors(c_pfas_diff3)
idxc3 <- rowSums(counts(c_pfas_diff3, normalized=TRUE) >= 5 ) >= 3
c_pfas_diff3 <- c_pfas_diff3[idxc3,]
dds_c_pfas3 <- DESeq(c_pfas_diff3)
res_c_pfas3 <- results(dds_c_pfas3)

#PFDA
c_pfas_diff4 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + mat_PFDA + cord_PFDA)
c_pfas_diff4 <- estimateSizeFactors(c_pfas_diff4)
idxc4 <- rowSums(counts(c_pfas_diff4, normalized=TRUE) >= 5 ) >= 3
c_pfas_diff4 <- c_pfas_diff4[idxc4,]
dds_c_pfas4 <- DESeq(c_pfas_diff4)
res_c_pfas4 <- results(dds_c_pfas4)

#PFUnDA
c_pfas_diff5 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + mat_PFUnDA + cord_PFUnDA)
c_pfas_diff5 <- estimateSizeFactors(c_pfas_diff5)
idxc5 <- rowSums(counts(c_pfas_diff5, normalized=TRUE) >= 5 ) >= 3
c_pfas_diff5 <- c_pfas_diff5[idxc5,]
dds_c_pfas5 <- DESeq(c_pfas_diff5)
res_c_pfas5 <- results(dds_c_pfas5)

#PFBS
c_pfas_diff6 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + mat_PFBS + cord_PFBS)
c_pfas_diff6 <- estimateSizeFactors(c_pfas_diff6)
idxc6 <- rowSums(counts(c_pfas_diff6, normalized=TRUE) >= 5 ) >= 3
c_pfas_diff6 <- c_pfas_diff6[idxc6,]
dds_c_pfas6 <- DESeq(c_pfas_diff6)
res_c_pfas6 <- results(dds_c_pfas6)

#PFHxS
c_pfas_diff7 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + mat_PFHxS + cord_PFHxS)
c_pfas_diff7 <- estimateSizeFactors(c_pfas_diff7)
idxc7 <- rowSums(counts(c_pfas_diff7, normalized=TRUE) >= 5 ) >= 3
c_pfas_diff7 <- c_pfas_diff7[idxc7,]
dds_c_pfas7 <- DESeq(c_pfas_diff7)
res_c_pfas7 <- results(dds_c_pfas7)


transcript <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.sex + condition.GA)
transcript <- estimateSizeFactors(transcript)
idx_trans <- rowSums(counts(transcript, normalized=TRUE) >= 5 ) >= 3
transcript <- transcript[idx_trans,]
diff_trans <- DESeq(transcript)
res_trans <- results(diff_trans)

save(res_c_pfas2, res_c_pfas3, res_c_pfas4, res_c_pfas5, res_c_pfas6, res_c_pfas7, res_c_pfas1, res_trans,
     file = "/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/pfas_results_taking_mat_control.RData")
