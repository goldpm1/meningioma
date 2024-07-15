#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long GENOME_FOLDER:,FASTQ_PATH_1:,FASTQ_PATH_2:,BAM_PATH:,OUTPUT_SORT_BAM:, -- "$@")
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
        --GENOME_FOLDER)
            GENOME_FOLDER=$2
        shift 2 ;;
        --FASTQ_PATH_1)
            FASTQ_PATH_1=$2
        shift 2 ;;
        --FASTQ_PATH_2)
            FASTQ_PATH_2=$2
        shift 2 ;;
        --BAM_PATH)
            BAM_PATH=$2
        shift 2 ;;
        --OUTPUT_SORT_BAM)
            OUTPUT_SORT_BAM=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e "bismark --genome ${GENOME_FOLDER} -1 ${FASTQ_PATH_1} -2 ${FASTQ_PATH_2}  -o ${BAM_PATH} --bowtie2 -p 8 --score_min L,0,-0.2 "



# bismark --genome ${GENOME_FOLDER} \
# -1 ${FASTQ_PATH_1} \
# -2 ${FASTQ_PATH_2}  \
# -o ${BAM_PATH} \
# --bowtie2 \
# -p 4 \
# --score_min L,0,-0.2  



OUTPUT_BAM=$(find "${BAM_PATH}" -name "*.bam" | head -n 1)
echo -e "samtools sort  ${OUTPUT_BAM} -o ${OUTPUT_SORT_BAM}"

samtools sort  ${OUTPUT_BAM} -o ${OUTPUT_SORT_BAM}
samtools index ${OUTPUT_SORT_BAM}