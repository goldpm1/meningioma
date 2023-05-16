#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long  Sample_ID:,CASE_BAM_PATH:,CONTROL_BAM_PATH:,OUTPUT_VCF_PATH:,REF:, -- "$@")
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
        --CASE_BAM_PATH)
            CASE_BAM_PATH=$2
        shift 2 ;;
        --CONTROL_BAM_PATH)
            CONTROL_BAM_PATH=$2
        shift 2 ;;
        --OUTPUT_VCF_PATH)
            OUTPUT_VCF_PATH=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo /home/goldpm1/tools/octopus/bin/octopus -R ${REF} -P 2 -I ${CONTROL_BAM_PATH}" "${CASE_BAM_PATH}" " -C cancer -N ${Sample_ID}"_Blood"  -o ${OUTPUT_VCF_PATH}

source ~/.bash_profile

INPUT_BAM=${CONTROL_BAM_PATH}" "${CASE_BAM_PATH} 
echo ${INPUT_BAM}

#1. Octopus tumor calling (n < 10)
/home/goldpm1/tools/octopus/bin/octopus -R ${REF} -I ${INPUT_BAM} -C cancer -N ${Sample_ID}"_Blood"  -o ${OUTPUT_VCF_PATH} \ 
--threads 16 --mask-tails 5  --mask-soft-clipped-bases --min-read-length 80 --ignore-unmapped-contigs \
--forest "/home/goldpm1/tools/octopus/resources/forests/germline.v0.7.4.forest.gz"  \
--annotations DP AD AF SB GQ BQ STRL STRP \
--somatics-only \
--filter-expression "QUAL < 20 | ADP < 10 | STRL < 20"