#!/usr/bin/env Rscript

# Librerías necesarias
library(dplyr)
library(data.table)
library(table1)

# Capturar argumentos desde la línea de comandos
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Uso: script.R <input_df> <input_genes> <output_dir>")
}

# Asignar argumentos a variables
input_df <- args[1]       # Path al archivo de variantes
input_genes <- args[2]    # Path al archivo de genes objetivo
output_dir <- args[3]     # Directorio de salida

# Leer archivos
df <- fread(input_df, header = TRUE)
genes <- fread(input_genes, header = TRUE)

# Renombrar columnas de df
colnames(df) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER",
                  "ALLELE","EFFECT","IMPACT","GENE","GENEID","FEATURE","FEATUREID","HGVS_C","HGVS_P","CDNA_POS",
                  "CADD_phred","Polyphen2_HVAR_pred","SIFT_pred","MutationTaster_pred","MutationAssessor_pred",
                  "ESP6500_EA_AF","ExAC_AF","ExAC_NFE_AF","ExAC_NFE_AC","ExAC_Adj_AF","1000Gp3_AF",
                  "gnomAD_exomes_NFE_AC","gnomAD_exomes_NFE_AF","gnomAD_genomes_NFE_AF","gnomAD_genomes_AF","gnomAD_exomes_AF",
                  "LRT_score","REVEL_score","clinvar_id","clinvar_clnsig","clinvar_trait",
                  "clinvar_review","clinvar_hgvs","clinvar_MedGen_id","clinvar_OMIM_id","clinvar_Orphanet_id","SpliceAI")

# Unir los genes con el dataframe de variantes
genes1 <- genes %>% left_join(df, by = c("Gene_name" = "GENE"), multiple = "all")

# Advertencia sobre relaciones many-to-many
message("¡Atención! Posible relación many-to-many detectada en el join.")

# Generar tabla de resumen
tab <- table1::table1(~ EFFECT | Gene_name, data = genes1)

# Definir paths de salida
output_table <- file.path(output_dir, "candidate-GENES.table")
output_rdata <- file.path(output_dir, "candidate-GENES.Rdata")

# Guardar archivos de salida
write.table(tab, output_table, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
save(genes1, file = output_rdata)

message("Proceso completado. Archivos guardados en: ", output_dir)
