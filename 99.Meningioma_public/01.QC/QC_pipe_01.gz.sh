#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,FASTQ_PATH_1:,FASTQ_PATH_2:,FASTQ_RAW_PATH_1:,FASTQ_RAW_PATH_2: -- "$@")
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
        --FASTQ_RAW_PATH_1)
            FASTQ_RAW_PATH_1=$2
        shift 2 ;;
        --FASTQ_RAW_PATH_2)
            FASTQ_RAW_PATH_2=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


gzip -f ${FASTQ_RAW_PATH_1} 
gzip -f ${FASTQ_RAW_PATH_2} 

mv ${FASTQ_RAW_PATH_1} ${FASTQ_PATH_1}
mv ${FASTQ_RAW_PATH_2} ${FASTQ_PATH_2}