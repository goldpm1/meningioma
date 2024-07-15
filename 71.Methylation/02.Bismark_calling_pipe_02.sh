#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long DSS_PATH:,CALL_DIR:, -- "$@")
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
        --DSS_PATH)
            DSS_PATH=$2
        shift 2 ;;
        --CALL_DIR)
            CALL_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


CALL_PATH_GZ=$(find "${CALL_DIR}" -name "*.cov.gz" | head -n 1)
gunzip ${CALL_PATH_GZ}

CALL_PATH=$(find "${CALL_DIR}" -name "*.cov*" | head -n 1)

echo -e "python3 02.Bismark_calling_pipe_02.py --CALL_PATH ${CALL_PATH} --DSS_PATH ${DSS_PATH}"
python3 02.Bismark_calling_pipe_02.py --CALL_PATH ${CALL_PATH} --DSS_PATH ${DSS_PATH}