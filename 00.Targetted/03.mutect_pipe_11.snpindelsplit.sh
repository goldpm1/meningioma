#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long INPUT_VCF_GZ:,OUTPUT_SNP_VCF:,OUTPUT_INDEL_VCF:, -- "$@")
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
        --INPUT_VCF_GZ)
            INPUT_VCF_GZ=$2
        shift 2 ;;
        --OUTPUT_SNP_VCF)
            OUTPUT_SNP_VCF=$2
        shift 2 ;;
        --OUTPUT_INDEL_VCF)
            OUTPUT_INDEL_VCF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



java -Xmx8g -jar "/opt/Yonsei/Picard/2.25.1/picard.jar" SplitVcfs \
      -I ${INPUT_VCF_GZ} \
      -SNP_OUTPUT ${OUTPUT_SNP_VCF} \
      -INDEL_OUTPUT ${OUTPUT_INDEL_VCF} \
      -STRICT false
