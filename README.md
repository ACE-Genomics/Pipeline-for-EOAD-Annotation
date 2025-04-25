# Variant Annotation Pipeline

Este script automatiza el proceso de anotación de variantes a partir de archivos VCF utilizando SnpEff y dbNSFP, además de filtrar genes de interés y generar análisis de portadores.

## 🛠️ Requisitos Cluster
- Java (>= 11)
- Singularity
- SnpEff + SnpSift
- R y los paquetes necesarios (para los scripts `.R`)
- Plink + Plink2

## 📦 Archivos requeridos
- Archivo VCF a anotar
- Lista de genes de interés
- Imagen Singularity de snpEff
- Scripts `.R`: `Filter_EffectVar.R`, `Create_ADAD_summ.R`, `Summary_varandcarriers.R`

## 🚀 Uso

sbatch EOAD_Annotation_JA.sh \
--inputvcf ../data/wes_joint_chr.snps.indels.g_recalibrated.vcf.gz \
--phenofile /nas/Genomica/02-Projects/2023_EOAD-clinical_VF-R2/202411_EOAD-AES.pheno \
--tag EOAD_TS \
--genes 202410_target_gene_list.txt \
--workdir /nas/Genomica/02-Projects/2025_EOAD-clinical_AVS/2023_EOAD-clinical-TESTS/prueba_automatized/
