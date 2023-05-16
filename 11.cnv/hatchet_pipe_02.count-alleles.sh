#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#export GRB_LICENSE_FILE=/opt/Yonsei/Gurobi/9.1.1/linux64/gurobi_nodes.lic

if ! options=$(getopt -o h --long TUMORS:,SAMPLES:,NORMAL:,SNPS:,REF:,OUTPUTNORMAL:,OUTPUTTUMORS:,OUTPUTSNPS:,SAMTOOLS:,BCFTOOLS:,MINCOV:,MAXCOV:,PROCESSES:,READQUALITY:,BASEQUALITY:, -- "$@")
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
        --SNPS)
            SNPS=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --OUTPUTNORMAL)
            OUTPUTNORMAL=$2
        shift 2 ;;
        --OUTPUTTUMORS)
            OUTPUTTUMORS=$2
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
        --READQUALITY)
            READQUALITY=$2
        shift 2 ;;
        --BASEQUALITY)
            BASEQUALITY=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


date

TUMORS_NEW=$(echo ${TUMORS} | sed 's/[{}]//g' | sed 's/,/ /g')
SAMPLES_NEW=$(echo ${SAMPLES} | sed 's/[{}]//g' | sed 's/,/ /g')
SNPS_NEW=${SNPS}"/*.vcf.gz"
# echo ${TUMORS_NEW}
# echo ${SAMPLES_NEW}
# echo ${SNPS_NEW}


echo "python3 -m hatchet count-alleles \
    --tumors ${TUMORS_NEW} --samples ${SAMPLES_NEW} --normal ${NORMAL} --snps ${SNPS_NEW} --reference ${REF} \
    --outputnormal ${OUTPUTNORMAL}  --outputtumors ${OUTPUTTUMORS} --outputsnps ${OUTPUTSNPS} \
    --samtools ${SAMTOOLS} --bcftools ${BCFTOOLS} --mincov ${MINCOV} --maxcov ${MAXCOV}  --processes ${PROCESSES} --readquality ${READQUALITY} --basequality ${BASEQUALITY}"

python3 -m hatchet count-alleles \
    --tumors ${TUMORS_NEW} --samples ${SAMPLES_NEW} --normal ${NORMAL} --snps ${SNPS_NEW} --reference ${REF} \
    --outputnormal ${OUTPUTNORMAL}  --outputtumors ${OUTPUTTUMORS} --outputsnps ${OUTPUTSNPS} \
    --samtools ${SAMTOOLS} --bcftools ${BCFTOOLS} --mincov ${MINCOV} --maxcov ${MAXCOV}  --processes ${PROCESSES}  --readquality ${READQUALITY} --basequality ${BASEQUALITY}
