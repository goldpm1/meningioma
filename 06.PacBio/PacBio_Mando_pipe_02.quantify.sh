#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long FASTQ_PATH:,GTF_FILE:,REF:,OUTPUT_DIR:, -- "$@")
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
        --REF)
            REF=$2
        shift 2 ;;
        --GTF_FILE)
            GTF_FILE=$2
        shift 2 ;;
        --FASTQ_PATH)
            FASTQ_PATH=$2
        shift 2 ;;
        --OUTPUT_DIR)
            OUTPUT_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e "python3 /home/goldpm1/tools/Mandalorion/Mando.py -p ${OUTPUT_DIR} -g ${GTF_FILE} -G ${REF} -f ${FASTQ_PATH}"
python3 /home/goldpm1/tools/Mandalorion/Mando.py -p ${OUTPUT_DIR} -g ${GTF_FILE} -G ${REF} -f ${FASTQ_PATH}