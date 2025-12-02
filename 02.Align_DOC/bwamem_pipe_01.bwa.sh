#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long FASTQ_PATH_1:,FASTQ_PATH_2:,PRE_BAM_PATH:,REF:,TMP_PATH:, -- "$@")
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
        --FASTQ_PATH_1)
            FASTQ_PATH_1=$2
        shift 2 ;;
        --FASTQ_PATH_2)
            FASTQ_PATH_2=$2
        shift 2 ;;
        --PRE_BAM_PATH)
            PRE_BAM_PATH=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --TMP_PATH)
            TMP_PATH=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e ${REF}"\n"${FASTQ_PATH_1}"\n"${FASTQ_PATH_2}

#BWA-mem / Picard
bwa mem -t 10 -M $REF ${FASTQ_PATH_1} ${FASTQ_PATH_2} | java -Xmx48g -jar /opt/Yonsei/Picard/2.26.4/picard.jar SortSam \
    SO=coordinate \
    INPUT=/dev/stdin \
    OUTPUT=${PRE_BAM_PATH} \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true \
    MAX_RECORDS_IN_RAM=1000000 \
    TMP_DIR=$TMP_PATH

date
echo "bwamem done"
