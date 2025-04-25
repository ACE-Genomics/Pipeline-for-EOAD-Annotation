# Variant Annotation Pipeline

Este script automatiza el proceso de anotación de variantes raras a partir de archivos VCF utilizando SnpEff y dbNSFP, además de filtrar genes de interés y generar análisis de portadores.

## 🛠️ Requisitos Cluster
- Java (>= 11)
- Singularity
- SnpEff + SnpSift
- R y los paquetes necesarios (para los scripts `.R`)
- Plink + Plink2

## 📦 Archivos requeridos
- Archivo VCF a anotar. El archivo debe usar el formato comprimido (.vcf.gz)
- Lista de genes de interés. Este archivo debe seguir el siguiente formato:
Gene_group	Gene_name	Canonical_Transcript
0	APP	ENST00000346798.8
0	PSEN1	ENST00000324501.10
0	PSEN2	ENST00000366783.8
1	GRN	ENST00000053867.8
1	MAPT	ENST00000262410.10
1	CHMP2B	ENST00000263780.9
1	FUS	ENST00000254108.12
- Imagen Singularity de snpEff
- Scripts `.R`: `Filter_EffectVar.R`, `Create_ADAD_summ.R`, `Summary_varandcarriers.R`

## 🚀 Uso

sbatch EOAD_Annotation_JA.sh \
--inputvcf /ruta/a/archivo.vcf.gz \
--phenofile /ruta/a/archivo.pheno \
--tag nombre_del_tag \
--genes /ruta/a/lista_de_genes.txt \
--workdir /ruta/a/directorio_de_trabajo
