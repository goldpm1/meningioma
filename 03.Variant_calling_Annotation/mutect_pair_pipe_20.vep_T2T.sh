#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long INPUT_VCF:,OUTPUT_VCF:, -- "$@")
then
    echo "ERROR: invalid options"
    exit 1
fi

eval set -- $options

while true; do
    case "$1" in
        -h|--help)
            echo "Usage"
        shift ;;
        --INPUT_VCF)
            INPUT_VCF=$2
        shift 2 ;;
        --OUTPUT_VCF)
            OUTPUT_VCF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

conda deactivate

# Using VEP cache with (recommended) FASTA sequence (the most efficient way)
/opt/Yonsei/ensembl-vep/101.0/vep \
    -i ${INPUT_VCF} \
    -o ${OUTPUT_VCF} \
    --fasta "/home/goldpm1/reference/chmT2T/Homo_sapiens-GCA_009914755.4-softmasked.fa.gz" \
    --format "vcf" \
    --vcf \
    --dir_cache "/home/goldpm1/reference/chmT2T" \
    --species homo_sapiens_gca009914755v4 \
    --cache_version 107 \
    --custom /home/goldpm1/reference/chmT2T/clinvar_20240624_GCA_009914755.4.vcf.gz,Clinvar.vcf.gz,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN  \
    --plugin AlphaMissense,file=/home/goldpm1/reference/chmT2T/AlphaMissense_T2T.sorted.tsv.gz \
    --domains --symbol --canonical --protein --biotype --uniprot --variant_class \
    --offline  --no_stats --force_overwrite


bgzip -c -f ${OUTPUT_VCF} > ${OUTPUT_VCF}".gz"
tabix -p vcf ${OUTPUT_VCF}".gz"



# /opt/Yonsei/ensembl-vep/101.0/vep \
#     -i "/data/project/Meningioma/01.PacBio/WGRS/03.vcf/02.DeepVariant/02.PASS/231101_Tumor.T2T/231101_Tumor.T2T.chr22.FMC.HF.RMBLACK.vcf" \
#     -o "/data/project/Meningioma/01.PacBio/WGRS/03.vcf/02.DeepVariant/02.PASS/231101_Tumor.T2T/231101_Tumor.T2T.chr22.FMC.HF.RMBLACK.vep.vcf" \
#     --fasta "/home/goldpm1/reference/chmT2T/Homo_sapiens-GCA_009914755.4-softmasked.fa.gz" \
#     --format "vcf" \
#     --vcf \
#     --dir_cache "/home/goldpm1/reference/chmT2T" \
#     --species homo_sapiens_gca009914755v4 \
#     --cache_version 107 \
#     --custom /home/goldpm1/reference/chmT2T/clinvar_20240624_GCA_009914755.4.vcf.gz,clinvar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN  \
#     --plugin AlphaMissense,file=/home/goldpm1/reference/chmT2T/AlphaMissense_T2T.sorted.tsv.gz \
#     --domains --symbol --canonical --protein --biotype --uniprot --variant_class \
#     --offline --no_stats --force_overwrite
