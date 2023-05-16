#!/bin/bash
#$ -S /bin/bash
#$ -cwd
# Basic Argument

if ! options=$(getopt -o h --long ID:,CONTROL_BAM_PATH:,CASE_BAM_PATH:,OUTPUT_PREFIX:,dbSNP:,TARGETS:, -- "$@")
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
        --ID)
            ID=$2
        shift 2 ;;
        --CONTROL_BAM_PATH)
            CONTROL_BAM_PATH=$2
        shift 2 ;;
        --CASE_BAM_PATH)
            CASE_BAM_PATH=$2
        shift 2 ;;
        --OUTPUT_PREFIX)
            OUTPUT_PREFIX=$2
        shift 2 ;;
        --dbSNP)
            dbSNP=$2
        shift 2 ;;
        --TARGETS)
            TARGETS=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

module load cnv_facets
echo -e "cnv_facets.R -t ${CASE_BAM_PATH} -n ${CONTROL_BAM_PATH} -vcf ${dbSNP} --targets ${TARGETS} -o ${OUTPUT_PREFIX} -g hg38"
cnv_facets.R -t ${CASE_BAM_PATH} -n ${CONTROL_BAM_PATH} -vcf ${dbSNP} --targets ${TARGETS} -o ${OUTPUT_PREFIX} -g hg38
cnv_facets.R -p ${OUTPUT_PREFIX}".csv.gz" -o ${OUTPUT_PREFIX}