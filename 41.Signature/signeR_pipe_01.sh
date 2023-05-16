#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,SIGNER_INPUT_VCF_DIR:,SIGNER_RESULT_DIR:, -- "$@")
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
        --SIGNER_INPUT_VCF_DIR)
            SIGNER_INPUT_VCF_DIR=$2
        shift 2 ;;
        --SIGNER_RESULT_DIR)
            SIGNER_RESULT_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e "/opt/Yonsei/R/4.2.0/bin/Rscript signeR_pipe_01.R \
    --SIGNER_INPUT_VCF_DIR ${SIGNER_INPUT_VCF_DIR} \
    --SIGNER_RESULT_DIR ${SIGNER_RESULT_DIR}"

/opt/Yonsei/R/4.2.0/bin/Rscript signeR_pipe_01.execute.R \
    --SIGNER_INPUT_VCF_DIR ${SIGNER_INPUT_VCF_DIR} \
    --SIGNER_RESULT_DIR ${SIGNER_RESULT_DIR} \
    --SIGNER_RESULT_Phat_PATH ${SIGNER_RESULT_DIR}"/signatures_Phat.tsv"


source /home/goldpm1/.bashrc
conda activate cnvpytor
python3 signeR_pipe_01.cosmic_fit.py \
    --SIGNER_RESULT_Phat_PATH ${SIGNER_RESULT_DIR}"/signatures_Phat.tsv" \
    --SIGNER_RESULT_cosmic_fit_PATH ${SIGNER_RESULT_DIR}"/signatures_cosmic_fit.txt"
    