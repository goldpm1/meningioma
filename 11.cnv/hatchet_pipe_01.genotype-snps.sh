#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#export GRB_LICENSE_FILE=/opt/Yonsei/Gurobi/9.1.1/linux64/gurobi_nodes.lic

if ! options=$(getopt -o h --long NORMAL:,REF:,OUTPUTSNPS:,SAMTOOLS:,BCFTOOLS:,MINCOV:,MAXCOV:,PROCESSES:, -- "$@")
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
        --NORMAL)
            NORMAL=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --OUTPUTSNPS)
            OUTPUTSNPS=$2
        shift 2 ;;
        --SAMTOOLS)
            SAMTOOLS=$2
        shift 2 ;;
        --BCFTOOLS)
            BCFTOOLS=$2
        shift 2 ;;
        --MINCOV)
            MINCOV=$2
        shift 2 ;;
        --MAXCOV)
            MAXCOV=$2
        shift 2 ;;
        --PROCESSES)
            PROCESSES=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo "~/miniconda3/envs/cnvpytor/bin/python3 -m hatchet genotype-snps \
     --normal ${NORMAL} --reference ${REF} --outputsnps ${OUTPUTSNPS}  --samtools ${SAMTOOLS} --bcftools ${BCFTOOLS} --mincov ${MINCOV} --maxcov ${MAXCOV}  --processes ${PROCESSES}"


#echo ${NORMAL//,/ }

source activate cnvpytor

~/miniconda3/envs/cnvpytor/bin/python3 -m hatchet genotype-snps \
     --normal ${NORMAL} --reference ${REF} --outputsnps ${OUTPUTSNPS}  --samtools ${SAMTOOLS} --bcftools ${BCFTOOLS} --mincov ${MINCOV} --maxcov ${MAXCOV}  --processes ${PROCESSES} 
