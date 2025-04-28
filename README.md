# Variant Annotation Pipeline

Este script automatiza el proceso de anotación de variantes raras a partir de archivos VCF utilizando SnpEff y dbNSFP, además de filtrar genes de interés y generar análisis de portadores.

## 🛠️ Requisitos Cluster
- Java (>= 11)
- Singularity
- SnpEff + SnpSift
- R y los paquetes necesarios (para los scripts `.R`)
- Plink + Plink2

## 📦 Archivos requeridos
A) - **archivo.vcf.gz**: Archivo VCF a anotar. El archivo debe usar el formato comprimido (.vcf.gz)
B) - **lista_de_genes.txt:** Lista de genes de interés. Este archivo debe seguir el siguiente formato:
  
| Gene_group | Gene_name | Canonical_Transcript    |
|------------|-----------|-------------------------|
| 0          | APP       | ENST00000346798.8       |
| 0          | PSEN1     | ENST00000324501.10      |
| 0          | PSEN2     | ENST00000366783.8       |
| 1          | GRN       | ENST00000053867.8       |
| 1          | MAPT      | ENST00000262410.10      |
| 1          | CHMP2B    | ENST00000263780.9       |
| 1          | FUS       | ENST00000254108.12      |
| 1          | TARDBP    | ENST00000040877.2       |
| 1          | TBK1      | ENST00000331710.10      |

El archivo de ejemplo con este formato tiene el nombre **202410_target_gene_list.txt**, pero puede ser modificado por el usuario para incluir nuevas versiones o genes adicionales según sea necesario.
C) - **20240429_LGM_list_grCh38.xlsx:** Archivo EXCEL de referencia con metainformación de variantes de AlzForum
  | n | Gene_name | input         | transcript                | strand | gDNA                 | Type | CHR_POS       | CHR_POS_REF_ALT       | cDNA      | protein | Source | Dx | AlzForum_category   |
|---|-----------|---------------|----------------------------|--------|----------------------|------|---------------|-----------------------|-----------|---------|--------|----|---------------------|
| 1 | PSEN2     | PSEN2:p.T122P  | NM_000447 (protein_coding) | +      | chr1:g.226885545A>C   | SNP  | chr1:226885545 | chr1:226885545:A:C     | c.364A>C  | p.T122P  | DIAN   | AD | likely pathogenic   |
| 2 | PSEN2     | PSEN2:p.N141Y  | NM_000447 (protein_coding) | +      | chr1:g.226885602A>T   | SNP  | chr1:226885602 | chr1:226885602:A:T     | c.421A>T  | p.N141Y  | DIAN   | AD | likely pathogenic   |
| 3 | PSEN2     | PSEN2:p.N141I  | NM_000447 (protein_coding) | +      | chr1:g.226885603A>T   | SNP  | chr1:226885603 | chr1:226885603:A:T     | c.422A>T  | p.N141I  | DIAN   | AD | pathogenic           |
  De momento la pipeline solo permite usar estas variantes de referencia para Alzforum (_update_ para futuras versiones)
D)

E) - Imagen Singularity de snpEff
F) - Scripts `.R`: `Filter_EffectVar.R`, `Create_ADAD_summ.R`, `Summary_varandcarriers.R`

## 🚀 Uso

sbatch EOAD_Annotation_JA.sh \
--inputvcf /ruta/a/archivo.vcf.gz \
--phenofile /ruta/a/archivo.pheno \
--tag nombre_del_tag \
--genes /ruta/a/lista_de_genes.txt \
--workdir /ruta/a/directorio_de_trabajo
