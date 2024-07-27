library(DESeq2)
library(RUVSeq)
library(MASS)
library(dplyr)
library(readxl)

load("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/salmon_output_concise.RData")
#gse <- summarizeToGene(se)
counts_se <- round(assay(se))
#counts <- round(assay(gse))
condition.GA <- coldata$condition.GA
dds <- DESeqDataSetFromMatrix(countData = counts_se, colData = coldata, design = ~ condition.GA)
dds <- DESeq(dds)
res <- results(dds)
sorted_res <- res[order(res$padj, decreasing = TRUE), ]
housekeeping_genes <- rownames(sorted_res)[1:1000]

set <- RUVg(counts_se, housekeeping_genes, k=1)
normalized_counts <- set$normalizedCounts
W_1 <- set$W
coldata_w1 <- cbind(coldata, D = W_1)

pfas_w1 <- read_xlsx("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/pfas_w1_adjusted.xlsx")
pfas_w1_cord <- pfas_w1[!is.na(pfas_w1$cord_PFBA), ]
pfas_w1_mat <- pfas_w1[!is.na(pfas_w1$mat_PFBA), ]
pfas_w1_both <- pfas_w1_mat[!is.na(pfas_w1_mat$cord_PFBA), ]

pfas_w1_cord_J <- paste0("J", pfas_w1_cord$condition.ID)
pfas_w1_mat_J <- paste0("J", pfas_w1_mat$condition.ID)
pfas_w1_both_J <- paste0("J", pfas_w1_both$condition.ID)
pfas_w1_both <- pfas_w1_both %>%
  mutate(r_PFBA = mat_PFBA / cord_PFBA)
pfas_w1_both <- pfas_w1_both %>%
  mutate(r_PFOA = mat_PFOA / cord_PFOA)
pfas_w1_both <- pfas_w1_both %>%
  mutate(r_PFNA = mat_PFNA / cord_PFNA)
pfas_w1_both <- pfas_w1_both %>%
  mutate(r_PFDA = mat_PFDA / cord_PFDA)
pfas_w1_both <- pfas_w1_both %>%
  mutate(r_PFUnDA = mat_PFUnDA / cord_PFUnDA)
pfas_w1_both <- pfas_w1_both %>%
  mutate(r_PFBS = mat_PFBS / cord_PFBS)
pfas_w1_both <- pfas_w1_both %>%
  mutate(r_PFHxS = mat_PFHxS / cord_PFHxS)


counts_filtered_cord <- counts_se[, colnames(counts_se) %in% pfas_w1_cord_J ]
counts_filtered_mat <- counts_se[, colnames(counts_se) %in% pfas_w1_mat_J ]
counts_filtered_both <- counts_se[, colnames(counts_se) %in% pfas_w1_both_J ]



#mat/cord PFBA
r_pfas_diff1 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + r_PFBA)
r_pfas_diff1 <- estimateSizeFactors(r_pfas_diff1)
idxr1 <- rowSums(counts(r_pfas_diff1, normalized=TRUE) >= 5 ) >= 3
r_pfas_diff1 <- r_pfas_diff1[idxr1,]
dds_r_pfas1 <- DESeq(r_pfas_diff1)
res_r_pfas1_isoform <- results(dds_r_pfas1)

#mat/cord PFOA
r_pfas_diff2 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + r_PFOA)
r_pfas_diff2 <- estimateSizeFactors(r_pfas_diff2)
idxr2 <- rowSums(counts(r_pfas_diff2, normalized=TRUE) >= 5 ) >= 3
r_pfas_diff2 <- r_pfas_diff2[idxr2,]
dds_r_pfas2 <- DESeq(r_pfas_diff2)
res_r_pfas2_isoform <- results(dds_r_pfas2)

#mat/cord PFNA
r_pfas_diff3 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + r_PFNA)
r_pfas_diff3 <- estimateSizeFactors(r_pfas_diff3)
idxr3 <- rowSums(counts(r_pfas_diff3, normalized=TRUE) >= 5 ) >= 3
r_pfas_diff3 <- r_pfas_diff3[idxr3,]
dds_r_pfas3 <- DESeq(r_pfas_diff3)
res_r_pfas3_isoform <- results(dds_r_pfas3)

#mat/cord PFDA
r_pfas_diff4 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + r_PFDA)
r_pfas_diff4 <- estimateSizeFactors(r_pfas_diff4)
idxr4 <- rowSums(counts(r_pfas_diff4, normalized=TRUE) >= 5 ) >= 3
r_pfas_diff4 <- r_pfas_diff4[idxr4,]
dds_r_pfas4 <- DESeq(r_pfas_diff4)
res_r_pfas4_isoform <- results(dds_r_pfas4)

#mat/cord PFUnDA
r_pfas_diff5 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + r_PFUnDA)
r_pfas_diff5 <- estimateSizeFactors(r_pfas_diff5)
idxr5 <- rowSums(counts(r_pfas_diff5, normalized=TRUE) >= 5 ) >= 3
r_pfas_diff5 <- r_pfas_diff5[idxr5,]
dds_r_pfas5 <- DESeq(r_pfas_diff5)
res_r_pfas5_isoform <- results(dds_r_pfas5)

#mat/cord PFBS
r_pfas_diff6 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + r_PFBS)
r_pfas_diff6 <- estimateSizeFactors(r_pfas_diff6)
idxr6 <- rowSums(counts(r_pfas_diff6, normalized=TRUE) >= 5 ) >= 3
r_pfas_diff6 <- r_pfas_diff6[idxr6,]
dds_r_pfas6 <- DESeq(r_pfas_diff6)
res_r_pfas6_isoform <- results(dds_r_pfas6)

#mat/cord PFHxS
r_pfas_diff7 <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + r_PFHxS)
r_pfas_diff7 <- estimateSizeFactors(r_pfas_diff7)
idxr7 <- rowSums(counts(r_pfas_diff7, normalized=TRUE) >= 5 ) >= 3
r_pfas_diff7 <- r_pfas_diff7[idxr7,]
dds_r_pfas7 <- DESeq(r_pfas_diff7)
res_r_pfas7_isoform <- results(dds_r_pfas7)

#cord PFBA
cord_pfas_diff1 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFBA)
cord_pfas_diff1 <- estimateSizeFactors(cord_pfas_diff1)
idx1 <- rowSums(counts(cord_pfas_diff1, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff1 <- cord_pfas_diff1[idx1,]
dds_cord_pfas1 <- DESeq(cord_pfas_diff1)
res_cord_pfas1_isoform <- results(dds_cord_pfas1)

#cord PFOA
cord_pfas_diff2 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFOA)
cord_pfas_diff2 <- estimateSizeFactors(cord_pfas_diff2)
idx2 <- rowSums(counts(cord_pfas_diff2, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff2 <- cord_pfas_diff2[idx2,]
dds_cord_pfas2 <- DESeq(cord_pfas_diff2)
res_cord_pfas2_isoform <- results(dds_cord_pfas2)

#cord PFNA
cord_pfas_diff3 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFNA)
cord_pfas_diff3 <- estimateSizeFactors(cord_pfas_diff3)
idx3 <- rowSums(counts(cord_pfas_diff3, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff3 <- cord_pfas_diff3[idx3,]
dds_cord_pfas3 <- DESeq(cord_pfas_diff3)
res_cord_pfas3_isoform <- results(dds_cord_pfas3)

#cord PFDA
cord_pfas_diff4 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFDA)
cord_pfas_diff4 <- estimateSizeFactors(cord_pfas_diff4)
idx4 <- rowSums(counts(cord_pfas_diff4, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff4 <- cord_pfas_diff4[idx4,]
dds_cord_pfas4 <- DESeq(cord_pfas_diff4)
res_cord_pfas4_isoform <- results(dds_cord_pfas4)

#cord PFUnDA
cord_pfas_diff5 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFUnDA)
cord_pfas_diff5 <- estimateSizeFactors(cord_pfas_diff5)
idx5 <- rowSums(counts(cord_pfas_diff5, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff5 <- cord_pfas_diff5[idx5,]
dds_cord_pfas5 <- DESeq(cord_pfas_diff5)
res_cord_pfas5_isoform <- results(dds_cord_pfas5)

#cord PFBS
cord_pfas_diff6 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFBS)
cord_pfas_diff6 <- estimateSizeFactors(cord_pfas_diff6)
idx6 <- rowSums(counts(cord_pfas_diff6, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff6 <- cord_pfas_diff6[idx6,]
dds_cord_pfas6 <- DESeq(cord_pfas_diff6)
res_cord_pfas6_isoform <- results(dds_cord_pfas6)

#cord PFHxS
cord_pfas_diff7 <- DESeqDataSetFromMatrix(countData = counts_filtered_cord, colData = pfas_w1_cord, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.GA + condition.sex + cord_PFHxS)
cord_pfas_diff7 <- estimateSizeFactors(cord_pfas_diff7)
idx7 <- rowSums(counts(cord_pfas_diff7, normalized=TRUE) >= 5 ) >= 3
cord_pfas_diff7 <- cord_pfas_diff7[idx7,]
dds_cord_pfas7 <- DESeq(cord_pfas_diff7)
res_cord_pfas7_isoform <- results(dds_cord_pfas7)

#mat PFBA
mat_pfas_diff1 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFBA)
mat_pfas_diff1 <- estimateSizeFactors(mat_pfas_diff1)
idx_m1 <- rowSums(counts(mat_pfas_diff1, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff1 <- mat_pfas_diff1[idx_m1,]
dds_mat_pfas1 <- DESeq(mat_pfas_diff1)
res_mat_pfas1_isoform <- results(dds_mat_pfas1)

#mat PFOA
mat_pfas_diff2 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFOA)
mat_pfas_diff2 <- estimateSizeFactors(mat_pfas_diff2)
idx_m2 <- rowSums(counts(mat_pfas_diff2, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff2 <- mat_pfas_diff2[idx_m2,]
dds_mat_pfas2 <- DESeq(mat_pfas_diff2)
res_mat_pfas2_isoform <- results(dds_mat_pfas2)

#mat PFNA
mat_pfas_diff3 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFNA)
mat_pfas_diff3 <- estimateSizeFactors(mat_pfas_diff3)
idx_m3 <- rowSums(counts(mat_pfas_diff3, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff3 <- mat_pfas_diff3[idx_m3,]
dds_mat_pfas3 <- DESeq(mat_pfas_diff3)
res_mat_pfas3_isoform <- results(dds_mat_pfas3)

#mat PFDA
mat_pfas_diff4 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFDA)
mat_pfas_diff4 <- estimateSizeFactors(mat_pfas_diff4)
idx_m4 <- rowSums(counts(mat_pfas_diff4, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff4 <- mat_pfas_diff4[idx_m4,]
dds_mat_pfas4 <- DESeq(mat_pfas_diff4)
res_mat_pfas4_isoform <- results(dds_mat_pfas4)

#mat PFUnDA
mat_pfas_diff5 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFUnDA)
mat_pfas_diff5 <- estimateSizeFactors(mat_pfas_diff5)
idx_m5 <- rowSums(counts(mat_pfas_diff5, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff5 <- mat_pfas_diff5[idx_m5,]
dds_mat_pfas5 <- DESeq(mat_pfas_diff5)
res_mat_pfas5_isoform <- results(dds_mat_pfas5)

#mat PFBS
mat_pfas_diff6 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFBS)
mat_pfas_diff6 <- estimateSizeFactors(mat_pfas_diff6)
idx_m6 <- rowSums(counts(mat_pfas_diff6, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff6 <- mat_pfas_diff6[idx_m6,]
dds_mat_pfas6 <- DESeq(mat_pfas_diff6)
res_mat_pfas6_isoform <- results(dds_mat_pfas6)

#mat PFHxS
mat_pfas_diff7 <- DESeqDataSetFromMatrix(countData = counts_filtered_mat, colData = pfas_w1_mat, design = ~ W_1 + condition.sex + condition.mother_ethnicity + condition.GROUP + condition.GA + mat_PFHxS)
mat_pfas_diff7 <- estimateSizeFactors(mat_pfas_diff7)
idx_m7 <- rowSums(counts(mat_pfas_diff7, normalized=TRUE) >= 5 ) >= 3
mat_pfas_diff7 <- mat_pfas_diff7[idx_m7,]
dds_mat_pfas7 <- DESeq(mat_pfas_diff7)
res_mat_pfas7_isoform <- results(dds_mat_pfas7)

transcript <- DESeqDataSetFromMatrix(countData = counts_filtered_both, colData = pfas_w1_both, design = ~ W_1 + condition.mother_ethnicity + condition.GROUP + condition.sex + condition.GA)
transcript <- estimateSizeFactors(transcript)
idx_trans <- rowSums(counts(transcript, normalized=TRUE) >= 5 ) >= 3
transcript <- transcript[idx_trans,]
diff_trans <- DESeq(transcript)
res_trans_isoform <- results(diff_trans)

save(res_cord_pfas1_isoform, res_cord_pfas2_isoform, res_cord_pfas3_isoform, res_cord_pfas4_isoform, res_cord_pfas5_isoform, res_cord_pfas6_isoform, res_cord_pfas7_isoform, res_mat_pfas1_isoform, res_mat_pfas2_isoform, res_mat_pfas3_isoform, res_mat_pfas4_isoform, res_mat_pfas5_isoform, res_mat_pfas6_isoform, res_mat_pfas7_isoform, 
     res_r_pfas1_isoform, res_r_pfas2_isoform, res_r_pfas3_isoform, res_r_pfas4_isoform, res_r_pfas5_isoform, res_r_pfas6_isoform, res_r_pfas7_isoform, res_trans_isoform,
     file = "/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/pfas_results_isoform.RData")



