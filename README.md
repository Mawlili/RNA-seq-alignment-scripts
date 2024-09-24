This repository contains the work flow of analyzing RNA expression data, specifically: star alignment, salmon quantification, tximeta import, and DESeq diffeq analysis.

Scripts were ran in the following order:

test_plink1.lsf -> Alignment_script_v2.lsf -> tximeta_import.lsf -> pfas_batch.qmd  -> PCA.R -> different PFAS association tests

Results and diagnostics of the analysis can be found in significant genes.zip. The most interpretable results were produced with the ones that took maternal level into control (DESeqDataSetFromMatrix(..., colData = ..., design = ~ W_1 + ... + mat_PFBA + cord_PFBA))
