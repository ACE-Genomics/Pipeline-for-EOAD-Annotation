#!/usr/bin/env Rscript

# Librerías necesarias
library(dplyr)
library(data.table)
library(table1)
library(stringr)
library(readxl)

# Capturar argumentos desde la línea de comandos
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Uso: script.R <input_rdata> <input_varlist> <output_dir>")
}

# Asignar argumentos a variables
input_rdata <- args[1]    # Path al archivo .Rdata de genes1
input_varlist <- args[2]  # Path al archivo Excel de variantes
output_dir <- args[3]     # Directorio de salida

# Cargar datos
load(input_rdata)
varlist <- read_excel(input_varlist, sheet = 1)

# Procesamiento de datos
genes1$ID <- paste(genes1$CHROM, genes1$POS, genes1$REF, genes1$ALT, sep = ":")
genes1$CHR_POS <- paste(genes1$CHROM, genes1$POS, sep = ":")
genes1$CHR_POS_REF_ALT <- paste(genes1$CHROM, genes1$POS, genes1$REF, genes1$ALLELE, sep = ":")

# Arreglar REVEL
genes1$REVEL_score <- substring(genes1$REVEL_score, 1, 5)

# Convertir columnas a numérico
columnas_a_convertir <- c("CADD_phred", "REVEL_score", "ESP6500_EA_AF", "1000Gp3_AF",
                          "ExAC_AF", "ExAC_NFE_AF", "ExAC_NFE_AC", "ExAC_Adj_AF",
                          "gnomAD_exomes_NFE_AC", "gnomAD_exomes_NFE_AF", "gnomAD_exomes_AF",
                          "gnomAD_genomes_NFE_AF", "gnomAD_genomes_AF")

#genes1[columnas_a_convertir] <- lapply(genes1[columnas_a_convertir], as.numeric)
genes1[, (columnas_a_convertir) := lapply(.SD, as.numeric), .SDcols = columnas_a_convertir]

# Clasificación de variantes
varlist_snp <- varlist[varlist$Type == "SNP",]

p1_unc <- varlist_snp[varlist_snp$AlzForum_category %in% c("uncertain significance", "Unclear Pathogenicity"),]
p1_pato <- varlist_snp[varlist_snp$AlzForum_category %in% c("pathogenic", "likely pathogenic"),]


genes1$P1 <- 0
genes1$P1[genes1$CHR_POS_REF_ALT %in% p1_unc$CHR_POS_REF_ALT] <- 1
genes1$P1[genes1$CHR_POS_REF_ALT %in% p1_pato$CHR_POS_REF_ALT] <- 2


# Puntuación CADD
genes1$P2 <- cut(genes1$CADD_phred, 
                 breaks = c(-Inf, 20, 25, 30, 35, 40, Inf), 
                 labels = c(0, 1, 2, 3, 4, 5), 
                 right = FALSE)
genes1$P2[is.na(genes1$P2)] <- 0
genes1$P2 <- as.numeric(as.character(genes1$P2))

#print(table(genes1$P2))

# Puntuación MAF
#print(genes1$gnomAD_genomes_AF)

#table(genes1$P3)

#añadimos que sea NA en cada condicion para que no chafe el valor anterior su ya cumplia una condicion anterior.
genes1$P3 <- NA 

# Definir los umbrales y sus puntajes asociados
thresholds <- c(0.00001, 0.0001, 0.001, 0.01)
scores <- c(4, 3, 2, 1)  # Puntajes asociados a cada umbral

# Obtener el valor máximo en cada fila, ignorando NAs
max_value <- pmax(genes1$gnomAD_genomes_AF, genes1$gnomAD_exomes_AF, na.rm = TRUE)

# Asignar puntajes según el valor máximo de la fila
for (i in seq_along(thresholds)) {
  threshold <- thresholds[i]
  
  genes1$P3[is.na(genes1$P3) & max_value < threshold] <- scores[i]
}

# Si ambos valores son NA, asignar el puntaje más alto (4)
genes1$P3[is.na(genes1$gnomAD_genomes_AF) & is.na(genes1$gnomAD_exomes_AF)] <- 3

# Si al menos un valor es mayor que 0.01, asignar 0
genes1$P3[is.na(genes1$P3) & (genes1$gnomAD_genomes_AF > 0.01 | genes1$gnomAD_exomes_AF > 0.01)] <- 0

#spliceai_score
# Esto asume que SpliceAI tiene siempre 11 campos separados por "|"
split_spliceAI <- do.call(rbind, strsplit(genes1$SpliceAI, "\\|"))

# Paso 2: Convertir los campos de interés a numéricos
# Índices: 3 = DS_AG, 4 = DS_AL, 5 = DS_DG, 6 = DS_DL
DS_AG <- as.numeric(split_spliceAI[, 3])
DS_AL <- as.numeric(split_spliceAI[, 4])
DS_DG <- as.numeric(split_spliceAI[, 5])
DS_DL <- as.numeric(split_spliceAI[, 6])

# Paso 3: Calcular P4: 1 si alguno > 0.2, 0 si no
# Inicializar todo en 0
genes1$P4 <- 0

# Ahora asignar 1 solo donde se cumple la condición
genes1$P4[DS_AG > 0.2 | DS_AL > 0.2 | DS_DG > 0.2 | DS_DL > 0.2] <- 1

#print(genes1)
#print(genes1$Gene_name)

# Clasificación final de variantes
genes1$Variant_class <- genes1$P1 + genes1$P2 + genes1$P3 + genes1$P4
test <- genes1[which(genes1$Gene_name == "APOE"),]
test2 <- table(test$Variant_class, test$Gene_name)

# Filtrado especial de TREM2
genes1_trem2 <- genes1[genes1$Gene_name == "TREM2",]
genes1_trem2 <- genes1_trem2[order(genes1_trem2$FEATUREID, decreasing = TRUE),]
genes1_trem2_uniq <- genes1_trem2[!duplicated(genes1_trem2$HGVS_C),]

genes2 <- rbind(genes1[genes1$Gene_name != "TREM2",], genes1_trem2_uniq)

# Tabla de variantes por gen
tab2 <- table(genes2$Variant_class, genes2$Gene_name)

# Filtrado de INDELs
ref_indels <- varlist[varlist$Type != "SNP",]
search_indels <- c("dup", "del", "ins", "-")
genes2_indels <- do.call(rbind, lapply(search_indels, function(x) {
  rbind(genes2[grep(x, genes2$HGVS_C),], genes2[grep(x, genes2$HGVS_P),])
}))
genes2_indels <- genes2_indels[!duplicated(genes2_indels),]
genes2_indelsA <- genes2_indels[genes2_indels$Gene_name %in% c("APP", "PSEN1", "MAPT"),]

# Lista de variantes
genes2_varlist <- genes2$CHR_POS_REF_ALT

# Definir paths de salida
output_table <- file.path(output_dir, "table_genes2.txt")
output_variants <- file.path(output_dir, "genes2_info_variants_and_score.txt")
output_varlist <- file.path(output_dir, "varlist_genes2.txt")
output_rdata <- file.path(output_dir, "candidate-GENES-filtered.Rdata")

# Guardar archivos
write.table(tab2, output_table, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(genes2, output_variants, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(genes2_varlist, output_varlist, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(genes1, genes1_trem2_uniq, genes2, genes2_varlist, genes2_indelsA, ref_indels, varlist, 
     file = output_rdata)

message("Proceso completado. Archivos guardados en: ", output_dir)
