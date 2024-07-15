#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long SIGPROFILER_INPUT_VCF_DIR:,SIGPROFILER_RESULT_COSMIC_DIR:,SIGPROFILER_RESULT_EXTRACT_DIR:, -- "$@")
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
        --SIGPROFILER_INPUT_VCF_DIR)
            SIGPROFILER_INPUT_VCF_DIR=$2
        shift 2 ;;
        --SIGPROFILER_RESULT_COSMIC_DIR)
            SIGPROFILER_RESULT_COSMIC_DIR=$2
        shift 2 ;;
        --SIGPROFILER_RESULT_EXTRACT_DIR)
            SIGPROFILER_RESULT_EXTRACT_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done




/opt/Yonsei/python/3.8.1/bin/python3 sigprofiler_pipe_01.cosmic_vcf.py \
    --SIGPROFILER_INPUT_VCF_DIR ${SIGPROFILER_INPUT_VCF_DIR} \
    --SIGPROFILER_RESULT_DIR ${SIGPROFILER_RESULT_COSMIC_DIR}