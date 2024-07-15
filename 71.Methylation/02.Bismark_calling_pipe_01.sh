#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long GENOME_FOLDER:,BAM_PATH:,CALL_DIR:, -- "$@")
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
        --GENOME_FOLDER)
            GENOME_FOLDER=$2
        shift 2 ;;
        --BAM_PATH)
            BAM_PATH=$2
        shift 2 ;;
        --CALL_DIR)
            CALL_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e "bismark_methylation_extractor   \
    --bedGraph \
    --buffer_size 150G \
    --ucsc \
    --cytosine_report \
    --genome_folder ${GENOME_FOLDER} \
    -o ${CALL_DIR} \
    ${BAM_PATH}"



bismark_methylation_extractor   \
    --bedGraph \
    --buffer_size 150G \
    --ucsc \
    --cytosine_report \
    --genome_folder ${GENOME_FOLDER} \
    -o ${CALL_DIR} \
    ${BAM_PATH}