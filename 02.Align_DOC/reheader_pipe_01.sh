#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,TISSUE:,FINAL_BAM_PATH:,FINAL_BAM_PATH_PREVIOUS:, -- "$@")
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
        --TISSUE)
            TISSUE=$2
        shift 2 ;;
        --FINAL_BAM_PATH)
            FINAL_BAM_PATH=$2
        shift 2 ;;
        --FINAL_BAM_PATH_PREVIOUS)
            FINAL_BAM_PATH_PREVIOUS=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


cp ${FINAL_BAM_PATH_PREVIOUS} ${FINAL_BAM_PATH}

echo -e "samtools view -h -o ${FINAL_BAM_PATH%bam}"sam" ${FINAL_BAM_PATH}"
samtools view -h -o ${FINAL_BAM_PATH%bam}"sam" ${FINAL_BAM_PATH}

echo -e "sed -i "s/${Sample_ID}_${TISSUE}/${Sample_ID}_2_${TISSUE}/g" "${FINAL_BAM_PATH%.bam}.sam""
sed -i "s/${Sample_ID}_${TISSUE}/${Sample_ID}_2_${TISSUE}/g" "${FINAL_BAM_PATH%.bam}.sam"

echo -e "samtools view -b -o ${FINAL_BAM_PATH%bam}"reheader.bam" ${FINAL_BAM_PATH%bam}"sam""
samtools view -b -o ${FINAL_BAM_PATH%bam}"reheader.bam" ${FINAL_BAM_PATH%bam}"sam"

echo -e "samtools index ${FINAL_BAM_PATH%bam}"reheader.bam""
samtools index ${FINAL_BAM_PATH%bam}"reheader.bam"