#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long SCRIPT_DIR:,OUTPUT_SBS96:,EXTRACTOR_DIR:, -- "$@")
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
        --SCRIPT_DIR)
            SCRIPT_DIR=$2
        shift 2 ;;
        --OUTPUT_SBS96)
            OUTPUT_SBS96=$2
        shift 2 ;;
        --EXTRACTOR_DIR)
            EXTRACTOR_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done




python3 ${SCRIPT_DIR}"/sigprofiler_pipe_14.DeNovoExtractor.py" \
    --OUTPUT_SBS96 ${OUTPUT_SBS96} \
    --EXTRACTOR_DIR ${EXTRACTOR_DIR}