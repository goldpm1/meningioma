#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long TMP_PATH:,REF:,INTERVAL:,INPUT_PATH:,OUTPUT_PATH:, -- "$@")
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
        --TMP_PATH)
            TMP_PATH=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --INTERVAL)
            INTERVAL=$2
        shift 2 ;;
        --INPUT_PATH)
            INPUT_PATH=$2
        shift 2 ;;
        --OUTPUT_PATH)
            OUTPUT_PATH=$2
        shift 2 ;;

        --)
            shift
            break
    esac
done



echo "1.WES : /opt/Yonsei/samtools/1.7/samtools mpileup -f ${REF} -l ${INTERVAL} -o ${OUTPUT_PATH} ${INPUT_PATH}"
/opt/Yonsei/samtools/1.7/samtools mpileup -f ${REF} -l ${INTERVAL} -o ${OUTPUT_PATH} ${INPUT_PATH}
