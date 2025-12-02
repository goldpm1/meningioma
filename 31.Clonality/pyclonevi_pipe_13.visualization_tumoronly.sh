#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long Sample_ID:,FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH:,FACETCNV_TO_PYCLONEVI_OUTPUT_PATH:,OUTPUT_VIS_FIG1:,OUTPUT_VIS_FIG2:,OUTPUT_VIS_DF:,OUTPUT_DIR1:,OUTPUT_DIR2:, -- "$@")
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
        --FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH)
            FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH=$2
        shift 2 ;;
        --FACETCNV_TO_PYCLONEVI_OUTPUT_PATH)
            FACETCNV_TO_PYCLONEVI_OUTPUT_PATH=$2
        shift 2 ;;
        --OUTPUT_VIS_FIG1)
            OUTPUT_VIS_FIG1=$2
        shift 2 ;;
        --OUTPUT_VIS_FIG2)
            OUTPUT_VIS_FIG2=$2
        shift 2 ;;
        --OUTPUT_VIS_DF)
            OUTPUT_VIS_DF=$2
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


echo -e "python3 pyclonevi_pipe_13.visualization_tumoronly.py \
        --Sample_ID ${Sample_ID} \
        --FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH} \
        --FACETCNV_TO_PYCLONEVI_OUTPUT_PATH ${FACETCNV_TO_PYCLONEVI_OUTPUT_PATH} \
        --OUTPUT_VIS_FIG1 ${OUTPUT_VIS_FIG1} \
        --OUTPUT_VIS_FIG2 ${OUTPUT_VIS_FIG2} \
        --OUTPUT_VIS_DF ${OUTPUT_VIS_DF} \
        --OUTPUT_DIR1 ${OUTPUT_DIR1} \
        --OUTPUT_DIR2 ${OUTPUT_DIR2}"


sleep 15s

python3 pyclonevi_pipe_13.visualization_tumoronly.py \
        --Sample_ID ${Sample_ID} \
        --FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH} \
        --FACETCNV_TO_PYCLONEVI_OUTPUT_PATH ${FACETCNV_TO_PYCLONEVI_OUTPUT_PATH} \
        --OUTPUT_VIS_FIG1 ${OUTPUT_VIS_FIG1} \
        --OUTPUT_VIS_FIG2 ${OUTPUT_VIS_FIG2} \
        --OUTPUT_VIS_DF ${OUTPUT_VIS_DF} \
        --OUTPUT_DIR1 ${OUTPUT_DIR1} \
        --OUTPUT_DIR2 ${OUTPUT_DIR2}