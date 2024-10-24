library(biomaRt)
bed_file <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_BRCA_BED/named_lifted_log2_input_TCGA_BRCA.bed")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_ids <- bed_file$gene_ID
gene_ids_clean <- sub("\\..*", "", gene_ids)
bed_file$gene_id_clean <- sub("\\..*", "", bed_file$gene_ID)


biomart_results <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = gene_ids_clean,
  mart = mart
)

merged_data <- merge(bed_file, biomart_results, by.x = "gene_id_clean", by.y = "ensembl_gene_id", all.x = TRUE)
protein_coding <- subset(merged_data, gene_biotype == "protein_coding")
non_protein_coding <- subset(merged_data, gene_biotype != "protein_coding")
save(protein_coding, file = "protein_coding.RData")
save(non_protein_coding, file = "non_protein_coding.RData")
protein_coding_filtered <- protein_coding[, !colnames(protein_coding) %in% "gene_id_clean"]
non_protein_coding_filtered <- non_protein_coding[, !colnames(non_protein_coding) %in% "gene_id_clean"]

write.table(protein_coding_filtered,
            file = "protein_coding.bed",
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write.table(non_protein_coding_filtered,
            file = "non_protein_coding.bed",
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
