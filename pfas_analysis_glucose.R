library(DESeq2)
library(MASS)
library(dplyr)
library(readxl)
library(apeglm)

load("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/salmon_data_concise.RData")
#filter out NA 
w1_glucose <- coldata_w1[!is.na(coldata_w1$condition.gdm_who_1999), ] 
#filter for common sample id
w1_glucose$condition.ID <- paste0("J", w1_glucose$condition.ID)
counts_glucose <- counts[, colnames(counts) %in% w1_glucose$condition.ID ]

#construct deseqdataset and perform deseq analysis

#PFBA
c_pfas_diff1 <- DESeqDataSetFromMatrix(countData = counts_glucose, colData = w1_glucose, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + condition.gdm_who_1999)
c_pfas_diff1 <- estimateSizeFactors(c_pfas_diff1)
idxc1 <- rowSums(counts(c_pfas_diff1, normalized=TRUE) >= 5 ) >= 3
c_pfas_diff1 <- c_pfas_diff1[idxc1,]
dds_c_pfas1 <- DESeq(c_pfas_diff1)
res_glucose <- lfcShrink(dds_c_pfas1, coef = "condition.gdm_who_1999_Yes_vs_No", type = "apeglm")

save(res_glucose,
     file = "/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/glucose_diffeq.RData")
