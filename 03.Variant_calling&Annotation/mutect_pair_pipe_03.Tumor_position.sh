#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,TUMOR_INTERVAL:,OUTPUT_FMC_HF_RMBLACK_PATH:,BCFTOOLS_MERGE_TXT:, -- "$@")
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
        --Sample_ID)
            Sample_ID=$2
        shift 2 ;;
        --TUMOR_INTERVAL)
            TUMOR_INTERVAL=$2
        shift 2 ;;
        --OUTPUT_FMC_HF_RMBLACK_PATH)
            OUTPUT_FMC_HF_RMBLACK_PATH=$2
        shift 2 ;;
        --BCFTOOLS_MERGE_TXT)
            BCFTOOLS_MERGE_TXT=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



echo -e "Sample_ID : "${Sample_ID}" OUTPUT_FMC_HF_RMBLACK_PATH :"${OUTPUT_FMC_HF_RMBLACK_PATH}" TUMOR_INTERVAL : "${TUMOR_INTERVAL}
echo -e "bcftools query -l ${OUTPUT_FMC_HF_RMBLACK_PATH} | bcftools view ${OUTPUT_FMC_HF_RMBLACK_PATH} | awk '{print $1"\t"$2-1"\t"$2"\t"$4"\t"$5}' > ${TUMOR_INTERVAL}"

# Use bcftools to extract the genomic coordinates from the VCF and convert it to BED
bcftools query -l ${OUTPUT_FMC_HF_RMBLACK_PATH} | bcftools view ${OUTPUT_FMC_HF_RMBLACK_PATH} | awk '{print $1"\t"$2-1"\t"$2"\t"$4"\t"$5}' > ${TUMOR_INTERVAL}".temp"
grep -v '#' ${TUMOR_INTERVAL}".temp" > ${TUMOR_INTERVAL}
rm -rf ${TUMOR_INTERVAL}".temp"

# bcftools merge txt 파일 생성
echo ${OUTPUT_FMC_HF_RMBLACK_PATH}".gz" > ${BCFTOOLS_MERGE_TXT}