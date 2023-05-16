#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#export GRB_LICENSE_FILE=/opt/Yonsei/Gurobi/9.1.1/linux64/gurobi_nodes.lic

if ! options=$(getopt -o h --long ARRAY:,BAFFILE:,TOTALCOUNTS:,REFVERSION:,OUTFILE:,MSR:,MTR:,PROCESSES:, -- "$@")
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
        --ARRAY)
            ARRAY=$2
        shift 2 ;;
        --BAFFILE)
            BAFFILE=$2
        shift 2 ;;
        --TOTALCOUNTS)
            TOTALCOUNTS=$2
        shift 2 ;;
        --REFVERSION)
            REFVERSION=$2
        shift 2 ;;
        --OUTFILE)
            OUTFILE=$2
        shift 2 ;;
        --MSR)
            MSR=$2
        shift 2 ;;
        --MTR)
            MTR=$2
        shift 2 ;;
        --PROCESSES)
            PROCESSES=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


date


echo -e "python3 -m hatchet combine-counts \
    --array ${ARRAY} --baffile ${BAFFILE} --totalcounts ${TOTALCOUNTS} --refversion ${REFVERSION}  \
    --outfile ${OUTFILE} \
    --msr ${MSR} --mtr ${MTR} --processes ${PROCESSES}"


python3 -m hatchet combine-counts \
    --array ${ARRAY} --baffile ${BAFFILE} --totalcounts ${TOTALCOUNTS} --refversion ${REFVERSION}  \
    --outfile ${OUTFILE} \
    --msr ${MSR} --mtr ${MTR} --processes ${PROCESSES} 
