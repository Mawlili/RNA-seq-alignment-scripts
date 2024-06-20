library(tximeta)
subdirs <- list.dirs("/rsrch5/home/epi/bhattacharya_lab/data/mapqtl/GUSTO/salmon_output", recursive = FALSE)
quant_files <- file.path(subdirs, "quant.sf")
sample_id <- paste0("J", 1001:1200)
coldata <- data.frame(quant_files, names = sample_id, stringsAsFactors=FALSE)
se <- tximeta(coldata)
save
