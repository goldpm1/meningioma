#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,FASTQ_PATH_1:,FASTQ_PATH_2:,FASTP_PATH_1:,FASTP_PATH_2: -- "$@")
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
        --FASTQ_PATH_1)
            FASTQ_PATH_1=$2
        shift 2 ;;
        --FASTQ_PATH_2)
            FASTQ_PATH_2=$2
        shift 2 ;;
        --FASTP_PATH_1)
            FASTP_PATH_1=$2
        shift 2 ;;
        --FASTP_PATH_2)
            FASTP_PATH_2=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e "SampleID: "$Sample_ID"\nFASTQ_PATH_1: "$FASTQ_PATH_1"\nFASTQ_PATH_2: "$FASTQ_PATH_2"\nFASTP_PATH_1: "$FASTP_PATH_1"\nFASTP_PATH_2: "$FASTP_PATH_2

fastp \
-i ${FASTQ_PATH_1} -I ${FASTQ_PATH_2}  \
-o ${FASTP_PATH_1} -O ${FASTP_PATH_2}  \
-p \
-w 5 \
-5 \
-3 \
--trim_poly_g \
--trim_poly_x \
--length_required 15 \
-y --detect_adapter_for_pe \
-h ${FASTP_PATH_1%/*}/${Sample_ID}'.html'



gzip -f ${FASTP_PATH_1}
gzip -f ${FASTP_PATH_2}