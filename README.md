# Variant Annotation Pipeline

This script automates the annotation process of rare variants from VCF files using SnpEff and dbNSFP, as well as filtering genes of interest and generating carrier analyses.

## 🛠️ Requirements Cluster
- Java (>= 11)
- Singularity
- SnpEff + SnpSift
- R y los paquetes necesarios (para los scripts `.R`)
- Plink + Plink2

## 📦 Input files required
A) - **archivo.vcf.gz:** VCF file to annotate. The file must be in compressed format (.vcf.gz)

B) - **lista_de_genes.txt**: List of genes of interest. This file must follow the following format:
  
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

The sample file following this format is named **202410_target_gene_list.txt**, but it can be modified by the user to include new versions or additional genes as needed.

C) - **20240429_LGM_list_grCh38.xlsx**: Reference EXCEL file containing variant meta-information from AlzForum.
  | n | Gene_name | input         | transcript                | strand | gDNA                 | Type | CHR_POS       | CHR_POS_REF_ALT       | cDNA      | protein | Source | Dx | AlzForum_category   |
|---|-----------|---------------|----------------------------|--------|----------------------|------|---------------|-----------------------|-----------|---------|--------|----|---------------------|
| 1 | PSEN2     | PSEN2:p.T122P  | NM_000447 (protein_coding) | +      | chr1:g.226885545A>C   | SNP  | chr1:226885545 | chr1:226885545:A:C     | c.364A>C  | p.T122P  | DIAN   | AD | likely pathogenic   |
| 2 | PSEN2     | PSEN2:p.N141Y  | NM_000447 (protein_coding) | +      | chr1:g.226885602A>T   | SNP  | chr1:226885602 | chr1:226885602:A:T     | c.421A>T  | p.N141Y  | DIAN   | AD | likely pathogenic   |
| 3 | PSEN2     | PSEN2:p.N141I  | NM_000447 (protein_coding) | +      | chr1:g.226885603A>T   | SNP  | chr1:226885603 | chr1:226885603:A:T     | c.422A>T  | p.N141I  | DIAN   | AD | pathogenic           |
 
For now, the pipeline only allows the use of these reference variants for AlzForum (update planned for future versions).

D) **archivo.pheno** with phenotypes for individuals included in the variant annotation analysis:
Tab-delimited file where 1 column MUST contain sample ID code with the column name "IID"

| FID         | IID         | PID | MID | Sex | Status | Age  | Center              | GWAs_ID |
|-------------|-------------|-----|-----|-----|--------|------|---------------------|---------|
| 22D28227751 | 22D28227751 | 0   | 0   | 2   | 2      | 71.6 | CUN_Pamplona_hist    | BCN396  |
| 22D28227797 | 22D28227797 | 0   | 0   | 2   | 2      | 64   | Clinic_hist          | 1729    |

PSP-sample.pheno serves as a sample file for this requirement

E) - Imagen Singularity de snpEff

F) - Scripts `.R`: `Filter_EffectVar.R`, `Create_ADAD_summ.R`, `Summary_varandcarriers.R`

## 🚀 Uso

sbatch EOAD_Annotation_JA.sh \
--inputvcf /ruta/a/archivo.vcf.gz \
--phenofile /ruta/a/archivo.pheno \
--tag nombre_del_tag \
--genes /ruta/a/lista_de_genes.txt \
--workdir /ruta/a/directorio_de_trabajo
