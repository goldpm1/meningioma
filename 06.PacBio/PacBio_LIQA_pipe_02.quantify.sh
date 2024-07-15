#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long BAM_PATH:,GTF_COMPATIBLE_FILE:,OUTPUT_DIR:, -- "$@")
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
        --BAM_PATH)
            BAM_PATH=$2
        shift 2 ;;
        --GTF_COMPATIBLE_FILE)
            GTF_COMPATIBLE_FILE=$2
        shift 2 ;;
        --OUTPUT_DIR)
            OUTPUT_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo -e "liqa -task quantify  -refgene ${GTF_COMPATIBLE_FILE}	-out ${OUTPUT_DIR} 	-max_distance 20 -f_weight 1 -bam ${BAM_PATH}"


liqa -task quantify \
		-refgene ${GTF_COMPATIBLE_FILE}\
		-out ${OUTPUT_DIR}\
		-max_distance 20\
		-f_weight 1\
		-bam ${BAM_PATH}