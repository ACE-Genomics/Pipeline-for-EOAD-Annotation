# Variant Annotation Pipeline

Este script automatiza el proceso de anotación de variantes raras a partir de archivos VCF utilizando SnpEff y dbNSFP, además de filtrar genes de interés y generar análisis de portadores.

## 🛠️ Requisitos Cluster
- Java (>= 11)
- Singularity
- SnpEff + SnpSift
- R y los paquetes necesarios (para los scripts `.R`)
- Plink + Plink2

## 📦 Archivos requeridos
- archivo.vcf.gz: Archivo VCF a anotar. El archivo debe usar el formato comprimido (.vcf.gz)
- lista_de_genes.txt: Lista de genes de interés. Este archivo debe seguir el siguiente formato:

El archivo de ejemplo con este formato tiene el nombre **202410_target_gene_list.txt**, pero puede ser modificado por el usuario para incluir nuevas versiones o genes adicionales según sea necesario.

- Imagen Singularity de snpEff
- Scripts `.R`: `Filter_EffectVar.R`, `Create_ADAD_summ.R`, `Summary_varandcarriers.R`

## 🚀 Uso

sbatch EOAD_Annotation_JA.sh \
--inputvcf /ruta/a/archivo.vcf.gz \
--phenofile /ruta/a/archivo.pheno \
--tag nombre_del_tag \
--genes /ruta/a/lista_de_genes.txt \
--workdir /ruta/a/directorio_de_trabajo
