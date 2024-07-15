#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long Sample_ID:,SEQUENZA_TO_PYCLONEVI_MATRIX_PATH:,SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH:,FACETCNV_TO_PYCLONEVI_MATRIX_PATH:,FACETCNV_TO_PYCLONEVI_OUTPUT_PATH:,OUTPUT_PATH_SHARED:,OUTPUT_DIR1:,OUTPUT_DIR2:, -- "$@")
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
        --SEQUENZA_TO_PYCLONEVI_MATRIX_PATH)
            SEQUENZA_TO_PYCLONEVI_MATRIX_PATH=$2
        shift 2 ;;
        --SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH)
            SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH=$2
        shift 2 ;;
        --FACETCNV_TO_PYCLONEVI_MATRIX_PATH)
            FACETCNV_TO_PYCLONEVI_MATRIX_PATH=$2
        shift 2 ;;
        --FACETCNV_TO_PYCLONEVI_OUTPUT_PATH)
            FACETCNV_TO_PYCLONEVI_OUTPUT_PATH=$2
        shift 2 ;;
        --OUTPUT_PATH_SHARED)
            OUTPUT_PATH_SHARED=$2
        shift 2 ;;
        --OUTPUT_DIR1)
            OUTPUT_DIR1=$2
        shift 2 ;;
        --OUTPUT_DIR2)
            OUTPUT_DIR2=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



source /home/goldpm1/.bashrc
conda activate cnvpytor

python3 pyclonevi_pipe_03.visualization.py \
        --Sample_ID ${Sample_ID} \
        --SEQUENZA_TO_PYCLONEVI_MATRIX_PATH ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} \
        --SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH ${SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH} \
        --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
        --FACETCNV_TO_PYCLONEVI_OUTPUT_PATH ${FACETCNV_TO_PYCLONEVI_OUTPUT_PATH} \
        --OUTPUT_PATH_SHARED ${OUTPUT_PATH_SHARED} \
        --OUTPUT_DIR1 ${OUTPUT_DIR1} \
        --OUTPUT_DIR2 ${OUTPUT_DIR2}