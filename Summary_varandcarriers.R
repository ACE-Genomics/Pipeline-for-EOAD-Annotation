## Cargar librerías necesarias
library(dplyr)
library(data.table)
library(table1)
library(stringr)
library(readxl)
library(openxlsx)

# Obtener los argumentos de la línea de comandos
args <- commandArgs(trailingOnly = TRUE)

# Asignar los argumentos a variables
pheno_file <- args[1]
genotypes_file <- args[2]
genes_file <- args[3]
gene_groups_file <- args[4]
output_prefix <- args[5]

## Cargar datos previos
load(genes_file)

genotypes <- fread(genotypes_file, h=T)
pheno <- fread(pheno_file, h=T)
gene_groups <- fread(gene_groups_file, h=T)


## Curación de datos
# Verificar si genes2 es un data.frame

genes2 <- as.data.frame(genes2)
genes2$CHR_POS_REF_ALT_REF <- paste0(genes2$CHR_POS_REF_ALT, "_", genes2$REF)
genes2 <- genes2[,-c(1)] %>% distinct()


# Invertir dosis en genotipos
genotypes <- as.data.frame(genotypes)
genotypes[, 7:ncol(genotypes)] <- apply(genotypes[, 7:ncol(genotypes)], 2, function(x) {
  ifelse(x == 2, 0, ifelse(x == 0, 2, x))
})
colnames(genotypes) <- sub("_.*", "", colnames(genotypes))
genotypes <- genotypes[, -c(1, 3:6)]  # Eliminar columnas innecesarias

##### INFO POR INDIVIDUO
temp <- genotypes[, -1]  
#print(temp)
temp <- data.frame(lapply(temp, function(x) {
  x[x == 2] <- 1
  return(as.numeric(x))
}))
temp$carrier_of_Nvar <- apply(temp, 1, function(x) sum(x, na.rm = TRUE))
rownames(temp) <- genotypes$IID  
#print(temp)
names(temp) <- gsub("\\.", ":", names(temp))

## Añadir la info al pheno file
pheno2 <- pheno





#########################
#Variantes con score 4-7#
#########################

list_score4to6 <- genes2[which(genes2$Variant_class < 7 & genes2$Variant_class > 3), c("CHR_POS_REF_ALT", "Gene_name")]

# print(genes2[which(genes2$Variant_class > 6),])
#print(list_score7to9)build in case there are no variants TOOO with condition
if (nrow(list_score4to6) > 0) {
temp4a6 <-temp[, colnames(temp) %in% list_score4to6$CHR_POS_REF_ALT, drop = FALSE]
temp4a6$carrier_of_Nvar_4to6 <- rowSums(temp4a6, na.rm = TRUE)
  temp4a6[is.na(temp4a6)] <- 0  # Convertir NA en 901

temp4a6 <- as.data.frame(lapply(temp4a6, function(col) trimws(as.character(col))), stringsAsFactors = FALSE)
colnames(temp4a6) <- gsub("\\.", ":", colnames(temp4a6))

rownames(temp4a6) <- rownames(temp)
  # Transformar presencia de SNPs en nombres de columna

temp4a6_pos <- as.data.frame(apply(temp4a6, 1, function(row) {
  colnames(temp4a6)[row == 1] <- colnames(temp4a6)[row == 1]
  row[row == 1] <- colnames(temp4a6)[row == 1]
  return(row)
}))

# Return to function
#temp7a9_pos <- as.data.frame(apply(temp7a9_pos, 1, function(row) {
 # colnames(temp7a9_pos)[row == 1] <- colnames(temp7a9_pos)[row == 1]
 # row[row == 1] <- colnames(temp7a9_pos)[row == 1]
 # return(row)
#}))
  # Trasponer el dataframe modificado para mantener la estructura original
  temp4a6_pos <- t(temp4a6_pos)
  # Convertir la matriz resultante de nuevo a dataframe
  temp4a6_pos <- as.data.frame(temp4a6_pos, stringsAsFactors = FALSE)
  # Asegurar que los rownames se mantienen y filtrar rowname columnas
  rownames(temp4a6_pos) <- rownames(temp4a6)
        temp4a6_pos <- temp4a6_pos[, -c(ncol(temp4a6_pos)), drop = FALSE]

temp4a6_pos <- as.data.frame(lapply(temp4a6_pos, function(col) trimws(as.character(col))), stringsAsFactors = FALSE)
colnames(temp4a6_pos) <- gsub("\\.", ":", colnames(temp4a6_pos))
rownames(temp4a6_pos) <- rownames(temp4a6)

  # Añadir la columna concatenada
  temp4a6_pos$SNP_ID <- apply(temp4a6_pos, 1, function(row) {
    paste(row[row != "0"], collapse = "_")
  })

# Reemplazo por nombres de genes
  sustituciones <- setNames(list_score4to6$Gene_name, list_score4to6$CHR_POS_REF_ALT)
  temp4a6_gene <- as.data.frame(lapply(temp4a6_pos, function(col) {
    sapply(col, function(val) {
      if (val %in% names(sustituciones)) return(sustituciones[[val]]) else return(val)
    })
  }), stringsAsFactors = FALSE)

  # Concatenar nombres de genes
temp4a6_gene <- as.data.frame(lapply(temp4a6_gene, function(col) trimws(as.character(col))), stringsAsFactors = FALSE)
colnames(temp4a6_gene) <- gsub("\\.", ":", colnames(temp4a6_pos))
rownames(temp4a6_gene) <- rownames(temp4a6_pos)

  temp4a6_gene$Gene_name <- apply(temp4a6_gene[, -ncol(temp4a6_gene), drop=FALSE], 1, function(row) {
    paste(row[row != "0"], collapse = "_")
  })

} else {
  # Si no hay variantes relevantes con score 7-9, crear un dataframe vacío o con valores predeterminados
  temp4a6 <- data.frame(matrix(NA, nrow = nrow(temp), ncol = 0))
  rownames(temp4a6) <- rownames(temp)
  temp4a6_gene <- data.frame(Gene_name = rep(NA, nrow(temp)))
  rownames(temp4a6_gene) <- rownames(temp)
  temp4a6_pos <- data.frame(SNP_ID = rep(NA, nrow(temp)))
  rownames(temp4a6_pos) <- rownames(temp)
}




#########################
#Variantes con score 7-9#
#########################

list_score7to9 <- genes2[which(genes2$Variant_class > 6), c("CHR_POS_REF_ALT", "Gene_name")]

# print(genes2[which(genes2$Variant_class > 6),])
#print(list_score7to9)build in case there are no variants TOOO with condition
if (nrow(list_score7to9) > 0) {
temp7a9 <- temp[, colnames(temp) %in% list_score7to9$CHR_POS_REF_ALT, drop = FALSE]
temp7a9$carrier_of_Nvar_7to9 <- rowSums(temp7a9, na.rm = TRUE)
  temp7a9[is.na(temp7a9)] <- 0  # Convertir NA en 901

temp7a9 <- as.data.frame(lapply(temp7a9, function(col) trimws(as.character(col))), stringsAsFactors = FALSE)
colnames(temp7a9) <- gsub("\\.", ":", colnames(temp7a9))

rownames(temp7a9) <- rownames(temp)
  # Transformar presencia de SNPs en nombres de columna

temp7a9_pos <- as.data.frame(apply(temp7a9, 1, function(row) {
  colnames(temp7a9)[row == 1] <- colnames(temp7a9)[row == 1]
  row[row == 1] <- colnames(temp7a9)[row == 1]
  return(row)
}))

# Return to function
#temp7a9_pos <- as.data.frame(apply(temp7a9_pos, 1, function(row) {
 # colnames(temp7a9_pos)[row == 1] <- colnames(temp7a9_pos)[row == 1]
 # row[row == 1] <- colnames(temp7a9_pos)[row == 1]
 # return(row)
#}))  
  # Trasponer el dataframe modificado para mantener la estructura original

  temp7a9_pos <- t(temp7a9_pos)
  # Convertir la matriz resultante de nuevo a dataframe
  temp7a9_pos <- as.data.frame(temp7a9_pos, stringsAsFactors = FALSE)

  # Asegurar que los rownames se mantienen y filtrar rowname columnas
  rownames(temp7a9_pos) <- rownames(temp7a9)
	temp7a9_pos <- temp7a9_pos[, -c(ncol(temp7a9_pos)), drop = FALSE]

temp7a9_pos <- as.data.frame(lapply(temp7a9_pos, function(col) trimws(as.character(col))), stringsAsFactors = FALSE)
colnames(temp7a9_pos) <- gsub("\\.", ":", colnames(temp7a9_pos))
rownames(temp7a9_pos) <- rownames(temp7a9)

  # Añadir la columna concatenada
  temp7a9_pos$SNP_ID <- apply(temp7a9_pos, 1, function(row) {
    paste(row[row != "0"], collapse = "_")
  })

# Reemplazo por nombres de genes
  sustituciones <- setNames(list_score7to9$Gene_name, list_score7to9$CHR_POS_REF_ALT)
  temp7a9_gene <- as.data.frame(lapply(temp7a9_pos, function(col) {
    sapply(col, function(val) {
      if (val %in% names(sustituciones)) return(sustituciones[[val]]) else return(val)
    })
  }), stringsAsFactors = FALSE)

  # Concatenar nombres de genes
temp7a9_gene <- as.data.frame(lapply(temp7a9_gene, function(col) trimws(as.character(col))), stringsAsFactors = FALSE)
colnames(temp7a9_gene) <- gsub("\\.", ":", colnames(temp7a9_pos))
rownames(temp7a9_gene) <- rownames(temp7a9_pos)

  temp7a9_gene$Gene_name <- apply(temp7a9_gene[, -ncol(temp7a9_gene), drop=FALSE], 1, function(row) {
    paste(row[row != "0"], collapse = "_")
  })

} else {
  # Si no hay variantes relevantes con score 7-9, crear un dataframe vacío o con valores predeterminados
  temp7a9 <- data.frame(matrix(NA, nrow = nrow(temp), ncol = 0))
  rownames(temp7a9) <- rownames(temp)
  temp7a9_gene <- data.frame(Gene_name = rep(NA, nrow(temp)))
  rownames(temp7a9_gene) <- rownames(temp)
  temp7a9_pos <- data.frame(SNP_ID = rep(NA, nrow(temp)))
  rownames(temp7a9_pos) <- rownames(temp)
}





# Variantes en Genes0
genes2 <- left_join(genes2, gene_groups,
by=c("Gene_name" = "Gene_name", "Canonical_Transcript" = "Canonical_Transcript"))

list_genes0 <- genes2[which(genes2$Gene_group==0),c("CHR_POS_REF_ALT","Gene_name")]

#print(temp)
temp_g0 <- temp[, which(colnames(temp) %in% c(list_genes0$CHR_POS_REF_ALT)), drop = FALSE]

temp_g0$carrier_of_Nvar_g0 <- apply(temp_g0, 1, function(x) sum(x, na.rm = TRUE))






####################################
# Variantes patogénicas en AlzForum#
###################################
alzf <- genes2[genes2$P1 == 2, c("CHR_POS_REF_ALT", "Gene_name")]

if (nrow(alzf) > 0) {  # Solo procesar si hay variantes
  temp_alzf <- temp[, colnames(temp) %in% alzf$CHR_POS_REF_ALT, drop = FALSE]
  temp_alzf[is.na(temp_alzf)] <- 0

  # Convertir valores en nombres de SNPs
  temp_alzf_pos <- as.data.frame(t(apply(temp_alzf, 1, function(row) {
    row[row == 1] <- names(row)[row == 1]
    return(row)
  })))
  #print(temp_alzf_pos)
  rownames(temp_alzf_pos) <- rownames(temp_alzf)

  # Concatenar SNP_ID y Gene_name para AlzForum
  temp_alzf_pos$SNP_ID <- apply(temp_alzf_pos, 1, function(row) {
    paste(row[row != "0"], collapse = "_")
  })

  sustituciones <- setNames(alzf$Gene_name, alzf$CHR_POS_REF_ALT)
  temp_alzf_gene <- as.data.frame(lapply(temp_alzf_pos, function(col) {
    sapply(col, function(val) {
      if (val %in% names(sustituciones)) return(sustituciones[[val]]) else return(val)
    })
  }), stringsAsFactors = FALSE)
  temp_alzf_gene$Gene_name <- apply(temp_alzf_gene[, -ncol(temp_alzf_gene), drop = FALSE], 1, function(row) {
    paste(row[row != "0"], collapse = "_")
  })

  rownames(temp_alzf_gene) <- rownames(temp_alzf_pos)
  temp_alzf_gene$IID <- rownames(temp_alzf_gene)
  temp_alzf_gene$SNP_ID <- temp_alzf_pos$SNP_ID

} else {
  # Si no hay variantes en AlzForum, devolver dataframes vacíos
  temp_alzf <- data.frame(matrix(NA, nrow = nrow(temp), ncol = 0))
  rownames(temp_alzf) <- rownames(temp)
  temp_alzf_gene <- data.frame(Gene_name = rep(NA, nrow(temp)))
  rownames(temp_alzf_gene) <- rownames(temp)
  temp_alzf_pos <- data.frame(SNP_ID = rep(NA, nrow(temp)))
  rownames(temp_alzf_pos) <- rownames(temp)

}


#print(temp4a7_pos)
#print(temp_alzf_gene)

############################### ## Add info to the Pheno_file
genotypes$carrier_of_Nvar <- temp$carrier_of_Nvar 
genotypes$carrier_of_Nvar_4to6 <- temp4a6$carrier_of_Nvar_4to6
genotypes$carrier_of_Nvar_7to9 <- temp7a9$carrier_of_Nvar_7to9 
genotypes$carrier_of_Nvar_g0 <- temp_g0$carrier_of_Nvar_g0 
genotypes$SNP_ID_score7a9 <- temp7a9_pos$SNP_ID 
genotypes$GeneName_score7a9 <- temp7a9_gene$Gene_name  
genotypes$SNP_ID_score4a6 <- temp4a6_pos$SNP_ID
genotypes$GeneName_score4a6 <- temp4a6_gene$Gene_name

#print(colnames(genotypes))
#LEFT_JOIN TOTAL PHENO
#if (!is.null(genotypes$SNP_ID_score7a9)) {
library(tidyr)
pheno3 <- pheno2 %>%
    inner_join(genotypes[,c("IID","carrier_of_Nvar","carrier_of_Nvar_g0","carrier_of_Nvar_7to9","GeneName_score7a9","SNP_ID_score7a9",
"carrier_of_Nvar_4to6", "GeneName_score4a6", "SNP_ID_score4a6")], by = c("IID" = "IID")) 

# Separar SNP_ID_score7a9 y GeneName_score7a9 si contienen múltiples valores concatenados con "_"
pheno3 <- pheno3 %>%
    separate_rows(SNP_ID_score7a9, GeneName_score7a9, sep = "_") %>%
    separate_rows(SNP_ID_score4a6, GeneName_score4a6, sep = "_")

if (!is.null(genotypes$SNP_ID_score7a9)) {
pheno3 <- pheno3 %>%
    left_join(genes2[, c("CHR_POS_REF_ALT", "P1", "Variant_class")], 
              by = c("SNP_ID_score7a9" = "CHR_POS_REF_ALT"))}

if (!is.null(genotypes$SNP_ID_score4a6)) {
pheno3 <- pheno3 %>%
    left_join(genes2[, c("CHR_POS_REF_ALT", "P1", "Variant_class")],
              by = c("SNP_ID_score4a6" = "CHR_POS_REF_ALT"))}

# Unificar columnas si existen ambas, pero con sufijos para diferenciarlas
if ("P1.x" %in% names(pheno3) & "P1.y" %in% names(pheno3)) {
    pheno3 <- pheno3 %>%
        rename(AlzForum_score7a9 = P1.x, AlzForum_score4a6 = P1.y)
} else if ("P1.x" %in% names(pheno3)) {
    pheno3 <- pheno3 %>% rename(AlzForum_score7a9 = P1.x)
} else if ("P1.y" %in% names(pheno3)) {
    pheno3 <- pheno3 %>% rename(AlzForum_score4a6 = P1.y)
}

# Unificar Variant_class.x y Variant_class.y si existen ambas
if ("Variant_class.x" %in% names(pheno3) & "Variant_class.y" %in% names(pheno3)) {
    pheno3 <- pheno3 %>%
        rename(Score_score7a9 = Variant_class.x, Score_score4a6 = Variant_class.y)
} else if ("Variant_class.x" %in% names(pheno3)) {
    pheno3 <- pheno3 %>% rename(Score_score7a9 = Variant_class.x)
} else if ("Variant_class.y" %in% names(pheno3)) {
    pheno3 <- pheno3 %>% rename(Score_score4a6 = Variant_class.y)
}

# Reasignar valores en la columna AlzForum_7a9 y AlzForum_4a7 si existen
if ("AlzForum_score7a9" %in% names(pheno3)) {
    pheno3 <- pheno3 %>%
        mutate(AlzForum_score7a9 = recode(AlzForum_score7a9, `1` = "Uncertain", `2` = "Pathogenic", `0` = "Missing or Conflicting"))
}

if ("AlzForum_score4a6" %in% names(pheno3)) {
    pheno3 <- pheno3 %>%
        mutate(AlzForum_score4a6 = recode(AlzForum_score4a6, `1` = "Uncertain", `2` = "Pathogenic", `0` = "Missing or Conflicting"))
}


#pivot longer en distintas filas
if ("AlzForum_score4a6" %in% names(pheno3) & "AlzForum_score7a9" %in% names(pheno3)) {
pheno3_processed <- pheno3 %>% pivot_longer(
    cols = c(SNP_ID_score7a9, GeneName_score7a9, AlzForum_score7a9, Score_score7a9,
             SNP_ID_score4a6, GeneName_score4a6, AlzForum_score4a6, Score_score4a6),
    names_to = c(".value", "Group"),
    names_sep = "_score"
  ) %>% distinct() %>%
filter(!is.na(SNP_ID))
  # Eliminar las filas donde SNP_ID sea NA (por si acaso)
} else if ("AlzForum_score4a6" %in% names(pheno3)) {
pheno3_processed <- pheno3 %>% pivot_longer(
    cols = c(SNP_ID_score4a6, GeneName_score4a6, AlzForum_score4a6, Score_score4a6),
    names_to = c(".value", "Group"),
    names_sep = "_score"
  ) %>% distinct() %>%
filter(!is.na(SNP_ID))

} else if ("AlzForum_score7a9" %in% names(pheno3)) {
pheno3_processed <- pheno3 %>% pivot_longer(
    cols = c(SNP_ID_score7a9, GeneName_score7a9, AlzForum_score7a9, Score_score7a9),
    names_to = c(".value", "Group"),
    names_sep = "_score"
  ) %>% distinct() %>%
filter(!is.na(SNP_ID))

} else {
pheno3_processed <- pheno3 %>% distinct()
}

#print(pheno3)

# Subset con variantes patogénicas

#print(temp_alzf_gene)
iid_alzf_car <- rownames(temp_alzf_gene[!is.na(temp_alzf_gene$Gene_name) & temp_alzf_gene$Gene_name != "", ])
#print(temp_alzf_gene$SNP_ID)

#print(temp_alzf_gene$IID)
#print(pheno3_processed)
cols_to_remove <- c("Group", "SNP_ID", "GeneName", "AlzForum", "Score")

if (!is.null(temp_alzf_gene$SNP_ID)) {
  pheno3_alzf <- pheno3_processed[, !colnames(pheno3_processed) %in% cols_to_remove] %>%
distinct() %>%
    left_join(temp_alzf_gene[, c("IID", "SNP_ID", "Gene_name")], by = c("IID" = "IID")) %>%
    filter(IID %in% iid_alzf_car) %>%
    left_join(genes2[, c("CHR_POS_REF_ALT", "Variant_class")], by = c("SNP_ID" = "CHR_POS_REF_ALT"))
} else {
  pheno3_alzf <- pheno3_processed[, !colnames(pheno3_processed) %in% cols_to_remove] %>%
distinct() %>%
    filter(IID %in% iid_alzf_car)}

#print(colnames(pheno3_alzf))

#ADD PATHOGENIC INFO
#names(pheno3_alzf) <- c("Centro", "Code_original", "Muestra", "Sex", "AAO", "APOE", "Nvariants", "Nvar_genes0", "SNP_ID_AlzF", "Gene_AlzF", "Score", "Nvar_score7to9", "Gene_score7a9", "SNP_ID_score7a9")
#print(pheno3_alzf)
## Add the pathogenic into the all phenotype (not only the score7to9) 

#pheno3$OnlyPathogenicity <- NA
#pheno3$SNP_ID_OnlyPatho <- NA
#pheno3$Score_OnlyPatho <- NA
#pheno3$GeneName_OnlyPatho <- NA


# Iterar sobre las filas de pheno3_alzf para asignar valores a la nueva columna
#for (sample_pathogenic in 1:nrow(pheno3_alzf)) {
#    x <- which(pheno3$Muestra_EOAD == pheno3_alzf$Muestra_EOAD[sample_pathogenic])
#    pheno3$SNP_ID_OnlyPatho[x] <- pheno3_alzf$SNP_ID[sample_pathogenic]
#pheno3$Score_OnlyPatho[x] <- genes2$Variant_class[sample_pathogenic]
#pheno3$GeneName_OnlyPatho[x] <- pheno3_alzf$Gene_name[sample_pathogenic]

    # Asignar valores a las nuevas columnas de acuerdo a pheno3_alzf
    # Asignar el valor 'Pathogenic' a la nueva columna Pathogenicity
#    pheno3$OnlyPathogenicity[x] <- "Pathogenic"
#}


#print(pheno3)
#pheno3 <- pheno3[,c("N","Centro","Code_original","Muestra_EOAD","Sex","AAO","APOE","Round","carrier_of_Nvar","carrier_of_Nvar_g0", "carrier_of_Nvar_7to9", "GeneName_score7a9", "SNP_ID_score7a9", "AlzForum", "Score")] 
#names(pheno3) <- c("N","Centro", "Code_original", "Muestra_EOAD", "Sex", "AAO", "APOE", "Round", "Nvariants", "Nvar_genes0", "Nvar_score7to9", "Gene_score7a9", "SNP_ID_score7a9","AlzForum", "Score")   

#print(pheno3)
#Obtain info per variant in another file
 # Traspongo la matrix de genotipos para ver cuantos carriers hay por cada variante 
genotypes2 <- as.data.frame(t(genotypes[,-1])) 
colnames(genotypes2) <- genotypes$IID 
genotypes2$CHR_POS_REF_ALT <- rownames(genotypes2) 
# Ahora uso el temp porque esta en 0/1 non-carrier/carrier 
temp2 <- temp %>% select(-carrier_of_Nvar) 
temp3 <- t(temp2) 
colnames(temp3) <- rownames(temp) 
temp3 <- as.data.frame(temp3) 
temp3$Nindiv_carriers <- apply(temp3, 1, sum) 
temp3$CHR_POS_REF_ALT <- rownames(temp3) 

genes3 <- genes2 %>% left_join(temp3[,c("CHR_POS_REF_ALT","Nindiv_carriers")], by="CHR_POS_REF_ALT") %>% left_join(genotypes2, by="CHR_POS_REF_ALT")  
# Exportar resultados

write.table(genes3, paste0(output_prefix, "FullTable_EOAD_variants.txt"), quote=FALSE, sep="\t", row.names = FALSE, col.names=TRUE) 
write.xlsx(genes3, paste0(output_prefix, "FullTable_EOAD_variants.xlsx"))  


write.table(pheno3_processed, paste0(output_prefix, "pheno_EOAD_carriers.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)
write.xlsx(pheno3_processed, paste0(output_prefix, "pheno_EOAD_carriers.xlsx"))

write.table(pheno3_alzf, paste0(output_prefix, "pheno_EOAD_carriersAlzF.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)
write.xlsx(pheno3_alzf, paste0(output_prefix, "pheno_EOAD_carriersAlzF.xlsx"))

save(genes3, genes2, genes2_varlist, genes2_indelsA,ref_indels, varlist,pheno3_alzf, pheno3, genotypes,list_genes0, list_score7to9, alzf,temp, 
temp_g0, temp7a9_pos, temp_alzf_gene, file=paste0(output_prefix, "filtered_genes_carriers_summ.Rdata")) 
