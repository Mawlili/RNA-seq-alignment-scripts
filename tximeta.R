#library
library(tximeta)
library(readxl)

#construct dataset with covariate and file path
cov <- read_excel("/rsrch5/home/epi/bhattacharya_lab/data/mapqtl/GUSTO/covariates/20220823-Full_200_RNAseq_covars_v2.xlsx")
subdirs <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/mapqtl/GUSTO/salmon_output", recursive = FALSE)
files <- file.path(subdirs, "quant.sf")
sample_id <- paste0("J", 1001:1200)
coldata <- data.frame(files, names = sample_id, condition = cov, stringsAsFactors=FALSE)

#import data to R with tximeta
se <- tximeta(coldata,countsFromAbundance = 'lengthScaledTPM',dropInfReps=TRUE)
suppressPackageStartupMessages(library(SummarizedExperiment))
colData(se)
save.image("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/out/salmon_data_concise.RData")
