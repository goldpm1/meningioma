#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long THREADS:,INTERVAL:,INPUT_BAM:,BAM_MOSDEPTH_PREFIX:, -- "$@")
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
        --THREADS)
            THREADS=$2
        shift 2 ;;
        --INTERVAL)
            INTERVAL=$2
        shift 2 ;;
        --INPUT_BAM)
            INPUT_BAM=$2
        shift 2 ;;
        --BAM_MOSDEPTH_PREFIX)
            BAM_MOSDEPTH_PREFIX=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


if [ ! -d ${BAM_MOSDEPTH_PREFIX%/*} ] ; then
    mkdir -p ${BAM_MOSDEPTH_PREFIX%/*}
fi


source /home/goldpm1/.bashrc
conda activate cnvpytor

echo -e "mosdepth -t ${THREADS}  --no-per-base  --by ${INTERVAL}  ${BAM_MOSDEPTH_PREFIX} ${INPUT_BAM}"
mosdepth -t ${THREADS}  --no-per-base  --by ${INTERVAL}  ${BAM_MOSDEPTH_PREFIX} ${INPUT_BAM}