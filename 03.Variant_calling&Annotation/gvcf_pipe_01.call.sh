#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if ! options=$(getopt -o h --long CASE_BAM_PATH:,INTERVAL:,REF:,OUTPUT_GVCF:, -- "$@")
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
        --INTERVAL)
            INTERVAL=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --OUTPUT_GVCF)
            OUTPUT_GVCF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

gatk --java-options "-Xmx48g" HaplotypeCaller  \
   -R ${REF} \
   -L ${INTERVAL} \
   -I ${CASE_BAM_PATH} \
   -O ${OUTPUT_GVCF} \
   -ERC GVCF \
   --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter

bgzip -c -f ${OUTPUT_GVCF} > ${OUTPUT_GVCF}".gz"
tabix -p vcf ${OUTPUT_GVCF}".gz"