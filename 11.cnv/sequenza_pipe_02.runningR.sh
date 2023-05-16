#!/bin/bash
#$ -S /bin/bash
#$ -cwd


if ! options=$(getopt -o h --long ID:,hg:,SEQUENZA_SMALL_SEQZ:,SEQUENZA_OUTPUT_DIR:, -- "$@")
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
        --ID)
            ID=$2
        shift 2 ;;
        --hg)
            hg=$2
        shift 2 ;;
        --SEQUENZA_SMALL_SEQZ)
            SEQUENZA_SMALL_SEQZ=$2
        shift 2 ;;
        --SEQUENZA_OUTPUT_DIR)
            SEQUENZA_OUTPUT_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e "Rscript sequenza_pipe_02.R --ID ${ID} --hg ${hg} --SEQUENZA_SMALL_SEQZ ${SEQUENZA_SMALL_SEQZ} --SEQUENZA_OUTPUT_DIR ${SEQUENZA_OUTPUT_DIR}"

Rscript "sequenza_pipe_02.R" --ID ${ID} --hg ${hg} --SEQUENZA_SMALL_SEQZ ${SEQUENZA_SMALL_SEQZ} --SEQUENZA_OUTPUT_DIR ${SEQUENZA_OUTPUT_DIR}
