#!/bin/bash
#$ -cwd
#$ -S /bin/bash

# INPUT_TSV="/data/project/Alzheimer/EM_cluster/old/pilot/04.EM_input/pyclone_vi/pyclone_vi_220610.tsv"
# OUTPUT_H5="/data/project/Alzheimer/EM_cluster/old/pilot/04.EM_input/pyclone_vi/pyclone_vi_220610.h5"
# OUTPUT_TSV="/data/project/Alzheimer/EM_cluster/old/pilot/04.EM_input/pyclone_vi/pyclone_vi_output_220610.tsv"

if ! options=$(getopt -o h --long INPUT_TSV:,OUTPUT_H5:,OUTPUT_TSV:, -- "$@")
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
        --INPUT_TSV)
            INPUT_TSV=$2
        shift 2 ;;
        --OUTPUT_H5)
            OUTPUT_H5=$2
        shift 2 ;;
        --OUTPUT_TSV)
            OUTPUT_TSV=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



source /home/goldpm1/.bashrc
conda activate cnvpytor
module switch HDF5/1.10.7 HDF5/1.12.0


rm -rf ${OUTPUT_H5} ${OUPTUT_TSV}

echo -e "pyclone-vi fit -i ${INPUT_TSV} -o ${OUTPUT_H5} -c 6 -d beta-binomial -r 20 "
pyclone-vi fit -i ${INPUT_TSV} -o ${OUTPUT_H5} -c 6 -d beta-binomial -r 20 
echo -e "\npyclone-vi fit done"
date



echo -e "pyclone-vi write-results-file -i ${OUTPUT_H5} -o ${OUTPUT_TSV}"
pyclone-vi write-results-file -i ${OUTPUT_H5} -o ${OUTPUT_TSV}
echo -e "\npyclone-vi write-results-file done"
date
