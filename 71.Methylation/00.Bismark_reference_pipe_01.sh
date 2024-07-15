#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long GENOME_FOLDER:, -- "$@")
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
        --)
            shift
            break
    esac
done


echo -e "bismark_genome_preparation ${GENOME_FOLDER}"

bismark_genome_preparation ${GENOME_FOLDER}