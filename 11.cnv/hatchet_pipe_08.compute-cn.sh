#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#export GRB_LICENSE_FILE=/opt/Yonsei/Gurobi/9.1.1/linux64/gurobi_nodes.lic

if ! options=$(getopt -o h --long INPUT_BINS_SEGMENTS:,RUNNINGDIR:, -- "$@")
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
        --INPUT_BINS_SEGMENTS)
            INPUT_BINS_SEGMENTS=$2
        shift 2 ;;
        --RUNNINGDIR)
            RUNNINGDIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


date

echo -e "python3 -m hatchet compute-cn -i ${INPUT_BINS_SEGMENTS} --runningdir ${RUNNINGDIR} --clones 2,5 --diploid"

python3 -m hatchet compute-cn -i ${INPUT_BINS_SEGMENTS} --runningdir ${RUNNINGDIR} --clones 2,5 --diploid