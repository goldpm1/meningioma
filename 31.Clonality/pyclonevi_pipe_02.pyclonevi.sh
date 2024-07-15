#!/bin/bash
#$ -cwd
#$ -S /bin/bash

#/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/pyclone_vi/

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
#module rm HDF5
module load HDF5/1.12.0
module switch HDF5/1.10.7 HDF5/1.12.0
#module switch HDF5/1.10.7 HDF5/1.14.2


rm -rf ${OUTPUT_H5} ${OUTPUT_TSV}

echo -e "pyclone-vi fit -i ${INPUT_TSV} -o ${OUTPUT_H5} -c 6 -d beta-binomial -r 20 "
pyclone-vi fit -i ${INPUT_TSV} -o ${OUTPUT_H5} -c 6 -d beta-binomial -r 20 
echo -e "\npyclone-vi fit done"
date



echo -e "pyclone-vi write-results-file -i ${OUTPUT_H5} -o ${OUTPUT_TSV}"
pyclone-vi write-results-file -i ${OUTPUT_H5} -o ${OUTPUT_TSV}
echo -e "\npyclone-vi write-results-file done"
date



python3 /data/project/Meningioma/script/31.Clonality/pyclonevi_pipe_02.sort.py \
    --INPUT_TSV ${OUTPUT_TSV} \
    --OUTPUT_TSV ${OUTPUT_TSV}