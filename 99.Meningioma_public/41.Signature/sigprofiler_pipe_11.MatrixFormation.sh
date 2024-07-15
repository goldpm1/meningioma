#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long SCRIPT_DIR:,RUN:,TISSUE:,SIGPROFILER_INPUT_VCF_DIR:,SIGPROFILER_INPUT_MATRIX_DIR:, -- "$@")
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
        --RUN)
            RUN=$2
        shift 2 ;;
        --TISSUE)
            TISSUE=$2
        shift 2 ;;
        --SIGPROFILER_INPUT_VCF_DIR)
            SIGPROFILER_INPUT_VCF_DIR=$2
        shift 2 ;;
        --SIGPROFILER_INPUT_MATRIX_DIR)
            SIGPROFILER_INPUT_MATRIX_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

python3 ${SCRIPT_DIR}"/sigprofiler_pipe_11.MatrixFormation.py" \
    --RUN ${RUN} \
    --TISSUE ${TISSUE} \
    --SIGPROFILER_INPUT_VCF_DIR ${SIGPROFILER_INPUT_VCF_DIR} \
    --SIGPROFILER_INPUT_MATRIX_DIR ${SIGPROFILER_INPUT_MATRIX_DIR} \

python3 ${SCRIPT_DIR}"/sigprofiler_pipe_12.MatrixGenerator.py" \
    --PROJECT "Meningioma" \
    --OUTPUT_DIR ${SIGPROFILER_INPUT_MATRIX_DIR} 