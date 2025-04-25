#!/bin/bash
#SBATCH --array=1
#SBATCH --job-name=AnnotationEOAD
#SBATCH --output=out/out.txt
#SBATCH --error=err/err2.txt
#SBATCH --time=23:59:00

####   Script to perform the annotation pipeline from VCF files with SnpEff and dbNSFP  ####
############################################################################################
# 2025/03/14
# AlejandroValenzuelaSeba


mkdir -p out err

# USAGE:
# sbatch run_annotation.sh <VCF_FILENAME> <TAG> <DATE_TARGET_GENES_FILE> <WORKDIR>

### ---- HELP & USAGE ---- ###
print_help() {
    echo ""
    echo "🔧 Usage:"
    echo "  sbatch EOAD_Annotation_JA.sh --inputvcf <VCF_FILE> --phenofile <PHENO_FILE> --tag <TAG_NAME> --genes <GENES_LIST> --workdir <WORKDIR>"
    echo ""
    echo "🧾 Parameters:"
    echo "  --inputvc     Path to input VCF file (can be .gz)"
    echo "  --phenofile   Path to phenotype file"
    echo "  --tag         Analysis tag (e.g., 2023_EOAD-clinical)"
    echo "  --genes       Path to genes of interest list"
    echo "  --workdir     Working directory for outputs"
    echo ""
    echo "ℹ️ Example:"
    echo "  sbatch EOAD_Annotation_JA.sh \\"
    echo "    --inputvcf /nas/.../wes_joint_chr.snps.g_recalibrated.vcf \\"
    echo "    --phenofile /nas/.../pheno.txt \\"
    echo "    --tag 2023_EOAD-clinical \\"
    echo "    --genes /nas/.../202410_target_gene_list.txt \\"
    echo "    --workdir /nas/.../2023_EOAD-clinical-TESTS/"
    echo ""
}

# If user asks for help
if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    print_help
    exit 0
fi

### ---- PARSE USER ARGUMENTS ---- ###
VCF=""
PHENO=""
TAG=""
GENESFILE=""
WORKDIR=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --inputvcf) VCF="$2"; shift ;;
        --phenofile) PHENO="$2"; shift ;;
        --tag) TAG="$2"; shift ;;
        --genes) GENESFILE="$2"; shift ;;
        --workdir) WORKDIR="$2"; shift ;;
        *) echo "❌ Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

### ---- VALIDATION ---- ###
if [[ -z "$VCF" || -z "$PHENO" || -z "$TAG" || -z "$GENESFILE" || -z "$WORKDIR" ]]; then
    echo "❗ Usage: sbatch EOAD_Annotation_JA.sh --inputvcf <VCF> --phenofile <PHENO> --tag <TAG> --genes <GENESFILE> --workdir <WORKDIR>"
    exit 1
fi

### ---- USER PATHS ---- ###
DATA="$(dirname "$VCF")/"
VCF=$(basename "$VCF" .gz)


#USER PATHS
#DATA="$1"
#VCF="$2"  # ej: wes_joint_chr.snps.indels.g_recalibrated.vcf
#PHENO="$3"
#TAG="$4"  # ej: 2023_EOAD-clinical
#GENESFILE="$5"
#WORKDIR="$6"  # ej: /nas/Genomica/02-Projects/2025_EOAD-clinical_AVS/2023_EOAD-clinical-TESTS

## INNER PATHS --> modify your paths here
#BIN
ORIG_DIR="$(pwd)"
SNPEFF="/nas/software/snpEff/snpEff.jar"
SNPSIFT="/nas/software/snpEff/SnpSift.jar"
PLINK="/nas/software/plink"
PLINK2="/nas/software/plink2"
JAVA="/usr/lib/jvm/java-23-openjdk-23.0.2.0.7-1.rolling.el8.x86_64/bin/java"
SINGULARITY="singularity run -B /nas:/nas -B /nas/software/snpEff/data/:/opt/snpEff/data/ /nas/Genomica/01-Data/00-Reference_files/03-Container/snpeff.simg"
DBNSFP="/nas/software/snpEff/data/dbNSFP4.7a.txt.gz"
ONEPERLINE="/nas/software/snpEff/scripts/vcfEffOnePerLine.pl"
GATK='singularity run --cleanenv -B /home:/home -B /nas:/nas -B /ruby:/ruby -B /greebo:/greebo /nas/usr/local/opt/gatk4.simg gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx16G"'
REF="/nas/Genomica/01-Data/00-Reference_files/02-GRCh38/00_Bundle/Homo_sapiens_assembly38.fasta"
BCFTOOLS="/nas/software/bcftools/bin/bcftools"
SLURMDIR="${WORKDIR}/slurm_jobs/"
mkdir -p ${SLURMDIR}

#VCFs ANBNOTATION PATHS
SPLITDIR="${WORKDIR}00-VCFsAnnot/"
mkdir -p ${SPLITDIR}
ANNOTDIR="${WORKDIR}01-snpEFF/"
echo ${ANNOTDIR}
mkdir -p ${ANNOTDIR}
#FILTERING VARIANTS PATHS
FILTERDIR="${WORKDIR}02-Candidate_genes/"
mkdir -p ${FILTERDIR}
CLASSIFICATION="${WORKDIR}03-Classif_ADAD/"
mkdir -p ${CLASSIFICATION}
#CARRIER IDENTIFICATION PATHS
GENODIR="${WORKDIR}04-variants-carriers/"
mkdir -p ${GENODIR}

########## You shouldn't modify anything from now on!! ##########

#SPLIT MULTIALLELIC AND ANNOTATE VARIANTS
echo -e "1 - Decompress and split multiallelic variants"

gunzip -c ${DATA}${VCF}.gz > ${SPLITDIR}${VCF}
${GATK} SelectVariants -V ${SPLITDIR}${VCF} -R ${REF} --exclude-filtered -O ${SPLITDIR}${VCF%.*}-PASS.vcf.gz

echo -e "2 - Run bcftools to normalize SNVs"
${BCFTOOLS} norm -Oz -m -any ${SPLITDIR}${VCF} --threads 10 -o ${SPLITDIR}${VCF%.*.*}-norm.vcf.gz &> ${SPLITDIR}OUT1
#cat OUT1 >> "COUNTS_Variants_counts_by_QC.txt"
${BCFTOOLS} norm -Oz -f ${REF} ${SPLITDIR}${VCF%.*.*}-norm.vcf.gz --threads 10 -o ${SPLITDIR}${VCF%.*.*}-norm-ref.vcf.gz &> ${SPLITDIR}OUT2
${BCFTOOLS} annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' ${SPLITDIR}${VCF%.*.*}-norm-ref.vcf.gz --threads 10 -o ${SPLITDIR}${VCF%.*.*}-norm-ref-annotate.vcf.gz
gzip -d ${SPLITDIR}${VCF%.*.*}-norm-ref-annotate.vcf.gz


## 01) Simplify the file
#CHROM POS ID REF ALT QUAL FILTER INFO
echo "3- Simplify the vcf file"

MAX_HEADER_LINES=1000
head -${MAX_HEADER_LINES} ${SPLITDIR}${VCF%.*.*}-norm-ref-annotate.vcf | grep "^##" > ${ANNOTDIR}${VCF%.*}.vcf;
grep -v "^##" ${SPLITDIR}${VCF%.*.*}-norm-ref-annotate.vcf | cut -f1-8  >> ${ANNOTDIR}${VCF%.*}.vcf

echo "4- Annotating with SnpEff"

## 2) Annotate with SnpEff
${JAVA} -Xmx15g -jar ${SNPEFF} -v GRCh38.105 ${ANNOTDIR}${VCF%.*}.vcf > ${ANNOTDIR}${VCF%.*}-snpEff.vcf
mv snpEff_genes.txt ${ANNOTDIR}/${VCF%.*}_snpEff_genes.txt; \
mv snpEff_summary.html ${ANNOTDIR}/${VCF%.*}_snpEff_genes.html

echo "4- Annotation with SnpEff completed"
echo "5- Using dbSNFP to annotate custom fields and filter information on variants"

## 3) Annotate with dbNSFP custom fields
${JAVA} -Xmx15g -jar ${SNPSIFT} dbnsfp -v -db ${DBNSFP} -f CADD_phred,Polyphen2_HVAR_pred,SIFT_pred,MutationTaster_pred,MutationAssessor_pred,ESP6500_EA_AF,ExAC_AF,ExAC_NFE_AF,ExAC_NFE_AC,ExAC_Adj_AF,1000Gp3_AF,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AF,gnomAD_genomes_NFE_AF,gnomAD_genomes_AF,gnomAD_exomes_AF,LRT_score,REVEL_score,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id ${ANNOTDIR}/${VCF%.*}-snpEff.vcf > ${ANNOTDIR}/${VCF%.*}-snpEff-dbnsfp.vcf

## 4) Split per line
cat ${ANNOTDIR}/${VCF%.*}-snpEff-dbnsfp.vcf | ${ONEPERLINE} | ${JAVA} -Xmx15g -jar ${SNPSIFT} extractFields -e "."  - CHROM POS ID REF ALT QUAL FILTER "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "dbNSFP_CADD_phred" "dbNSFP_Polyphen2_HVAR_pred" "dbNSFP_SIFT_pred" "dbNSFP_MutationTaster_pred" "dbNSFP_MutationAssessor_pred" "dbNSFP_ESP6500_EA_AF" "dbNSFP_ExAC_AF" "dbNSFP_ExAC_NFE_AF" "dbNSFP_ExAC_NFE_AC" "dbNSFP_ExAC_Adj_AF" "dbNSFP_1000Gp3_AF" "dbNSFP_gnomAD_exomes_NFE_AC" "dbNSFP_gnomAD_exomes_NFE_AF" "dbNSFP_gnomAD_genomes_NFE_AF" "dbNSFP_gnomAD_genomes_AF" "dbNSFP_gnomAD_exomes_AF" "dbNSFP_LRT_score" "dbNSFP_REVEL_score" "dbNSFP_clinvar_id" "dbNSFP_clinvar_clnsig" "dbNSFP_clinvar_trait" "dbNSFP_clinvar_review" "dbNSFP_clinvar_hgvs" "dbNSFP_clinvar_MedGen_id" "dbNSFP_clinvar_OMIM_id" "dbNSFP_clinvar_Orphanet_id" > ${ANNOTDIR}/${VCF%.*}-snpEff-dbnsfp-FIELDS.txt

echo "5- Annotation of fields completed"

echo "6- Removing flagged variants and extracting variants and genes of interest from table"
# 7) Remove flagged variants and extracting variants and genes of interest
# First get the header
head -1 ${ANNOTDIR}${VCF%.*}-snpEff-dbnsfp-FIELDS.txt > ${ANNOTDIR}${VCF%.*}-snpEff-dbnsfp-FIELDS.txt-header

# 8.a) Remove non-PASS variants --> only for vcfs that still have this coding
grep -w "PASS" ${ANNOTDIR}${VCF%.*}-snpEff-dbnsfp-FIELDS.txt > ${ANNOTDIR}${VCF%.*}-snpEff-dbnsfp-FIELDS.txt-PASS # este paso es esencial cuando trabajamos con los EOAD files que no les hicimos el VQSR o plink filtering standard

## 8.b) Extracting variants of interest ##
while IFS=$'\t' read -r col1 col2 col3 || [[ -n "$col1" ]]; do
    # Genera la orden para cada fila
   echo -e "$col1\t$col2\t$col3";\
   echo -e "grep -w \"$col2\" ${ANNOTDIR}${VCF%.*}-snpEff-dbnsfp-FIELDS.txt | grep -w \"$col3\" | grep -E \"HIGH|MODERATE\" >> ${FILTERDIR}/SUBSET_${col2}_${col3}" >> ${FILTERDIR}pre_$(date +%y%m%d)_gene_filter_parallel.txt;\
done < 202410_target_gene_list.txt

sed '1d' ${FILTERDIR}pre_$(date +%y%m%d)_gene_filter_parallel.txt > ${FILTERDIR}$(date +%y%m%d)_gene_filter_parallel.txt
rm ${FILTERDIR}pre_$(date +%y%m%d)_gene_filter_parallel.txt

## 9) Run the script from the same directory that we have been creating it
cp ${FILTERDIR}$(date +%y%m%d)_gene_filter_parallel.txt ${SLURMDIR}
cp 00_slurmize_v3.pl ${SLURMDIR}
cd ${SLURMDIR}
./00_slurmize_v3.pl $(date +%y%m%d)_gene_filter_parallel.txt
sleep 60
cd ${ORIG_DIR}

echo "6- Removal and extraction of variants completed"

echo "7- Join all variants from different genes"
## 10) put all the selected variants "together". Accounting for duplicates (Vicky extra code)
cat ${FILTERDIR}SUBSET_* | sort | uniq > ${FILTERDIR}${VCF%.*}-candidate-GENES.txt
grep -w "PASS" ${FILTERDIR}${VCF%.*}-candidate-GENES.txt > ${FILTERDIR}${VCF%.*}-candidate-GENES.txt-PASS

echo "7- Joined variants"
## 11) Script to filter variants by effect (Rscript)

pwd
echo "8- Filtering variants by effect and summary classfification following ADAD"
Rscript Filter_EffectVar.R ${FILTERDIR}${VCF%.*}-candidate-GENES.txt-PASS \
202410_target_gene_list.txt \
${CLASSIFICATION}

## 12) Create classification summary (Rscript)

Rscript Create_ADAD_summ.R ${CLASSIFICATION}candidate-GENES.Rdata \
20240429_LGM_list_grCh38.xlsx \
${CLASSIFICATION}

echo "8- Filtered and classfied with ADAD"
echo "9- Identifying carriers"

## 13) identify carriers
#First account for the possiubility of duplicates

cat ${FILTERDIR}SUBSET_* > ${FILTERDIR}SUBSET_ALL_GENES_INTEREST.txt
cut -f 3 ${FILTERDIR}SUBSET_ALL_GENES_INTEREST.txt | sort | uniq > ${FILTERDIR}SUBSET_ALL_GENES_INTEREST.variants

#Generate plink files
${PLINK} --vcf ${SPLITDIR}${VCF%.*.*}-norm-ref-annotate.vcf --double-id --make-bed --keep-allele-order --out ${GENODIR}EOAD-vars-carriers

# Let's calcualte a MAF of 1 AC per 1064 people --> 1 allele between 2128 alelles = 0.00005
${PLINK} --bfile ${GENODIR}EOAD-vars-carriers --extract ${FILTERDIR}SUBSET_ALL_GENES_INTEREST.txt --make-bed --keep-allele-order --maf 0.00005 --out ${GENODIR}EOAD-vars-carriers-maf-LIST2

# Get summary stats for these variants
${PLINK2} --bfile ${GENODIR}EOAD-vars-carriers-maf-LIST2 --freq --hardy --geno-counts --sample-counts --out ${GENODIR}EOAD-vars-STATS-LIST2

# Recode into raw file
${PLINK2} --bfile ${GENODIR}EOAD-vars-carriers-maf-LIST2 --recode A --out ${GENODIR}EOAD-vars-carriers-maf-LIST2

#Create Excel summary variants and carriers

Rscript Summary_varandcarriers.R ${PHENO} \
${GENODIR}EOAD-vars-carriers-maf-LIST2.raw \
${CLASSIFICATION}candidate-GENES-filtered.Rdata \
${GENODIR}${TAG}

echo "9- Summary created. Pipeline completed"

