#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long GTF_FILE:,GTF_COMPATIBLE_FILE:, -- "$@")
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
        --GTF_FILE)
            GTF_FILE=$2
        shift 2 ;;
        --GTF_COMPATIBLE_FILE)
            GTF_COMPATIBLE_FILE=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo -e "liqa -task refgene -format gtf -ref ${GTF_FILE} -out ${GTF_COMPATIBLE_FILE}"
liqa -task refgene -format gtf -ref ${GTF_FILE} -out ${GTF_COMPATIBLE_FILE}