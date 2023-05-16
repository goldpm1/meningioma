#!/bin/bash

if ! options=$(getopt -o h --long Sample_ID:,FASTQ_PATH_1:,FASTQ_PATH_2:,FASTQC_PATH:, -- "$@")
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
        --FASTQC_PATH)
            FASTQC_PATH=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e "Sample_ID: "$Sample_ID"\nFASTQ_PATH_1: "$FASTQ_PATH_1"\nFASTQ_PATH_2: "$FASTQ_PATH_2"\nFASQTC_PATH: "$FASTQC_PATH


fastqc -o ${FASTQC_PATH} \
            -f fastq ${FASTQ_PATH_1} ${FASTQ_PATH_2}