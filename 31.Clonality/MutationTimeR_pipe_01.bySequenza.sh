#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,TISSUE:,SEQUENZA_SEGMENT_PATH:,SEQUENZA_PLOIDY_PATH:,MUTATIONTIMER_INPUT_VCF_PATH:,MUTATIONTIMER_RESULT_DIR:, -- "$@")
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
        --SEQUENZA_SEGMENT_PATH)
            SEQUENZA_SEGMENT_PATH=$2
        shift 2 ;;
        --SEQUENZA_PLOIDY_PATH)
            SEQUENZA_PLOIDY_PATH=$2
        shift 2 ;;
        --MUTATIONTIMER_INPUT_VCF_PATH)
            MUTATIONTIMER_INPUT_VCF_PATH=$2
        shift 2 ;;
        --MUTATIONTIMER_RESULT_DIR)
            MUTATIONTIMER_RESULT_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e "/opt/Yonsei/R/4.2.0/bin/Rscript MutationTimeR_pipe_01.bySequenza.R \
    --Sample_ID ${Sample_ID} \
    --TISSUE ${TISSUE} \
    --SEQUENZA_SEGMENT_PATH ${SEQUENZA_SEGMENT_PATH} \
    --SEQUENZA_PLOIDY_PATH ${SEQUENZA_PLOIDY_PATH} \
    --MUTATIONTIMER_INPUT_VCF_PATH ${MUTATIONTIMER_INPUT_VCF_PATH} \
    --MUTATIONTIMER_RESULT_DIR ${MUTATIONTIMER_RESULT_DIR}
"


/opt/Yonsei/R/4.2.0/bin/Rscript MutationTimeR_pipe_01.bySequenza.R \
    --Sample_ID ${Sample_ID} \
    --TISSUE ${TISSUE} \
    --SEQUENZA_SEGMENT_PATH ${SEQUENZA_SEGMENT_PATH} \
    --SEQUENZA_PLOIDY_PATH ${SEQUENZA_PLOIDY_PATH} \
    --MUTATIONTIMER_INPUT_VCF_PATH ${MUTATIONTIMER_INPUT_VCF_PATH} \
    --MUTATIONTIMER_RESULT_DIR ${MUTATIONTIMER_RESULT_DIR}
    