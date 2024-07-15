#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if ! options=$(getopt -o h --long CASE_BAM_PATH:,WES_TUMOR_BED_PATH:,OUTPUT_BAMSNAP_DIR:,SAMPLE_ID:,BAM_DIR_LIST:,TITLE_LIST:,REF:,OUTPUT_SCCALLER:, -- "$@")
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
        --CASE_BAM_PATH)
            CASE_BAM_PATH=$2
        shift 2 ;;
        --WES_TUMOR_BED_PATH)
            WES_TUMOR_BED_PATH=$2
        shift 2 ;;
        --OUTPUT_BAMSNAP_DIR)
            OUTPUT_BAMSNAP_DIR=$2
        shift 2 ;;
        --SAMPLE_ID)
            SAMPLE_ID=$2
        shift 2 ;;
        --BAM_DIR_LIST)
            BAM_DIR_LIST=$2
        shift 2 ;;
        --TITLE_LIST)
            TITLE_LIST=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --OUTPUT_SCCALLER)
            OUTPUT_SCCALLER=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


#1. sccaller 돌리기
echo -e "python2 "/home/goldpm1/tools/SCcaller/sccaller_v2.0.0.py" --bam ${CASE_BAM_PATH} --fasta ${REF} --output ${OUTPUT_SCCALLER} --snp_type hsnp --snp_in ${WES_TUMOR_BED_PATH} --cpu_num 16 --engine samtools"
python2 "/home/goldpm1/tools/SCcaller/sccaller_v2.0.0.py" --bam ${CASE_BAM_PATH} --fasta ${REF} --output ${OUTPUT_SCCALLER} --snp_type hsnp --snp_in ${WES_TUMOR_BED_PATH} --cpu_num 16 --engine samtools


# 2. Bamsnap 찍기
source /home/goldpm1/.bashrc
conda activate cnvpytor