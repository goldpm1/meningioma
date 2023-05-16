#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#export GRB_LICENSE_FILE=/opt/Yonsei/Gurobi/9.1.1/linux64/gurobi_nodes.lic

if ! options=$(getopt -o h --long TUMORS:,SAMPLES:,NORMAL:,BAFFILE:,REFVERSION:,OUTDIR:,SAMTOOLS:,PROCESSES:,READQUALITY:, -- "$@")
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
        --TUMORS)
            TUMORS=$2
        shift 2 ;;
        --SAMPLES)
            SAMPLES=$2
        shift 2 ;;
        --NORMAL)
            NORMAL=$2
        shift 2 ;;
        --BAFFILE)
            BAFFILE=$2
        shift 2 ;;
        --REFVERSION)
            REFVERSION=$2
        shift 2 ;;
        --OUTDIR)
            OUTDIR=$2
        shift 2 ;;
        --SAMTOOLS)
            SAMTOOLS=$2
        shift 2 ;;
        --PROCESSES)
            PROCESSES=$2
        shift 2 ;;
        --READQUALITY)
            READQUALITY=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


date

TUMORS_NEW=$(echo ${TUMORS} | sed 's/[{}]//g' | sed 's/,/ /g')
SAMPLES_NEW=$(echo ${SAMPLES} | sed 's/[{}]//g' | sed 's/,/ /g')
# echo ${TUMORS_NEW}
# echo ${SAMPLES_NEW}


echo -e "python3 -m hatchet count-reads \
    --tumors ${TUMORS_NEW} --samples ${SAMPLES_NEW} --normal ${NORMAL} --refversion ${REFVERSION} --baffile ${BAFFILE} \
    --outdir ${OUTDIR} \
    --samtools ${SAMTOOLS} --processes ${PROCESSES}  --readquality ${READQUALITY} "

python3 -m hatchet count-reads \
    --tumors ${TUMORS_NEW} --samples ${SAMPLES_NEW} --normal ${NORMAL} --refversion ${REFVERSION} --baffile ${BAFFILE} \
    --outdir ${OUTDIR} \
    --samtools ${SAMTOOLS} --processes ${PROCESSES}  --readquality ${READQUALITY} 
