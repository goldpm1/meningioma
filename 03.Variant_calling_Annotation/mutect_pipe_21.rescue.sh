#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,TISSUE:,NUM:,MULTIPLE_VCF_PATH:,MULTIPLE_VCF_GZ_PATH:,INDIVIDUAL_VCF_PATH:,INDIVIDUAL_VCF_GZ_PATH:,INDIVIDUAL_RESCUED_VCF_PATH:,INDIVIDUAL_RESCUED_VCF_GZ_PATH:,INDIVIDUAL_UNIQUE_VCF_PATH:,INDIVIDUAL_UNIQUE_VCF_GZ_PATH:, -- "$@")
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
        --NUM)
            NUM=$2
        shift 2 ;;
        --MULTIPLE_VCF_PATH)
            MULTIPLE_VCF_PATH=$2
        shift 2 ;;
        --MULTIPLE_VCF_GZ_PATH)
            MULTIPLE_VCF_GZ_PATH=$2
        shift 2 ;;
        --INDIVIDUAL_VCF_PATH)
            INDIVIDUAL_VCF_PATH=$2
        shift 2 ;;
        --INDIVIDUAL_VCF_GZ_PATH)
            INDIVIDUAL_VCF_GZ_PATH=$2
        shift 2 ;;
        --INDIVIDUAL_RESCUED_VCF_PATH)
            INDIVIDUAL_RESCUED_VCF_PATH=$2
        shift 2 ;;
        --INDIVIDUAL_RESCUED_VCF_GZ_PATH)
            INDIVIDUAL_RESCUED_VCF_GZ_PATH=$2
        shift 2 ;;
        --INDIVIDUAL_UNIQUE_VCF_PATH)
            INDIVIDUAL_UNIQUE_VCF_PATH=$2
        shift 2 ;;
        --INDIVIDUAL_UNIQUE_VCF_GZ_PATH)
            INDIVIDUAL_UNIQUE_VCF_GZ_PATH=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


source /home/goldpm1/.bashrc
conda activate cnvpytor

echo -e "python3 mutect_pipe_21.rescue.py \
    --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} --NUM ${NUM} \
    --MULTIPLE_VCF_PATH ${MULTIPLE_VCF_PATH}  --MULTIPLE_VCF_GZ_PATH ${MULTIPLE_VCF_GZ_PATH} \
    --INDIVIDUAL_VCF_PATH ${INDIVIDUAL_VCF_PATH} --INDIVIDUAL_VCF_GZ_PATH ${INDIVIDUAL_VCF_GZ_PATH} \
    --INDIVIDUAL_RESCUED_VCF_PATH ${INDIVIDUAL_RESCUED_VCF_PATH} --INDIVIDUAL_RESCUED_VCF_GZ_PATH ${INDIVIDUAL_RESCUED_VCF_GZ_PATH} \
    --INDIVIDUAL_UNIQUE_VCF_PATH ${INDIVIDUAL_UNIQUE_VCF_PATH} --INDIVIDUAL_UNIQUE_VCF_GZ_PATH ${INDIVIDUAL_UNIQUE_VCF_GZ_PATH}"


python3 mutect_pipe_21.rescue.py \
    --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} --NUM ${NUM} \
    --MULTIPLE_VCF_PATH ${MULTIPLE_VCF_PATH}  --MULTIPLE_VCF_GZ_PATH ${MULTIPLE_VCF_GZ_PATH} \
    --INDIVIDUAL_VCF_PATH ${INDIVIDUAL_VCF_PATH} --INDIVIDUAL_VCF_GZ_PATH ${INDIVIDUAL_VCF_GZ_PATH} \
    --INDIVIDUAL_RESCUED_VCF_PATH ${INDIVIDUAL_RESCUED_VCF_PATH} --INDIVIDUAL_RESCUED_VCF_GZ_PATH ${INDIVIDUAL_RESCUED_VCF_GZ_PATH} \
    --INDIVIDUAL_UNIQUE_VCF_PATH ${INDIVIDUAL_UNIQUE_VCF_PATH} --INDIVIDUAL_UNIQUE_VCF_GZ_PATH ${INDIVIDUAL_UNIQUE_VCF_GZ_PATH}
