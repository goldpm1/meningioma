#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long INPUT_VCF:,OUTPUT_MAF:,SELECTED_DB:, -- "$@")
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
        --INPUT_VCF)
            INPUT_VCF=$2
        shift 2 ;;
        --OUTPUT_MAF)
            OUTPUT_MAF=$2
        shift 2 ;;
        --SELECTED_DB)
            SELECTED_DB=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


CURRENT_PATH=`pwd -P`
echo -e "python3 mutect_pipe_23.manualmaf.py  --INPUTVCF "${INPUT_VCF}" --OUTPUT_MAF "${OUTPUT_MAF}" --SELECTED_DB "${SELECTED_DB}
python3 "mutect_pipe_23.manualmaf.py" --INPUTVCF ${INPUT_VCF} --OUTPUT_MAF ${OUTPUT_MAF} --SELECTED_DB ${SELECTED_DB}