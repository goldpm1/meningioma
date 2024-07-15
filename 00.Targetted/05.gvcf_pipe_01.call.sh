#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if ! options=$(getopt -o h --long INPUT_BAM_PATH:,INTERVAL:,REF:,OUTPUT_GVCF_PATH:, -- "$@")
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
        --INPUT_BAM_PATH)
            INPUT_BAM_PATH=$2
        shift 2 ;;
        --INTERVAL)
            INTERVAL=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --OUTPUT_GVCF_PATH)
            OUTPUT_GVCF_PATH=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo -e "gatk --java-options "-Xmx48g" HaplotypeCaller  \
   -R ${REF} \
   -L ${INTERVAL} \
   -I ${INPUT_BAM_PATH} \
   -O ${OUTPUT_GVCF_PATH} \
   -ERC GVCF \
   --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter"

gatk --java-options "-Xmx48g" HaplotypeCaller  \
   -R ${REF} \
   -L ${INTERVAL} \
   -I ${INPUT_BAM_PATH} \
   -O ${OUTPUT_GVCF_PATH} \
   -ERC GVCF \
   --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter

bgzip -c -f ${OUTPUT_GVCF_PATH} > ${OUTPUT_GVCF_PATH}".gz"
tabix -f -p vcf ${OUTPUT_GVCF_PATH}".gz"