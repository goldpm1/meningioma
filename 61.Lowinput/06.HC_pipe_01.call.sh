#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if ! options=$(getopt -o h --long BAM_PATH:,INTERVAL:,REF:,dbSNP:,OUTPUT_HC_GVCF:,OUTPUT_HC_VCF:, -- "$@")
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
        --BAM_PATH)
            BAM_PATH=$2
        shift 2 ;;
        --INTERVAL)
            INTERVAL=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --dbSNP)
            dbSNP=$2
        shift 2 ;;
        --OUTPUT_HC_GVCF)
            OUTPUT_HC_GVCF=$2
        shift 2 ;;
        --OUTPUT_HC_VCF)
            OUTPUT_HC_VCF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


gatk --java-options "-Xmx48g" HaplotypeCaller  \
   -R ${REF} \
   -L ${INTERVAL} \
   -I ${BAM_PATH} \
   -D ${dbSNP} \
   --dont-use-soft-clipped-bases true \
    -ERC GVCF \
    --standard-min-confidence-threshold-for-calling 10 \
   -O ${OUTPUT_HC_GVCF} 

# 추후 MutationTimeR를 위해서 FORMAT 을 바꿔주기
sed 's/ID=AD,Number=R/ID=AD,Number=2/' ${OUTPUT_HC_GVCF} > ${OUTPUT_HC_GVCF}".temp"
mv ${OUTPUT_HC_GVCF}".temp" ${OUTPUT_HC_GVCF}
rm -rf ${OUTPUT_HC_GVCF}".temp"

bgzip -c -f ${OUTPUT_HC_GVCF} > ${OUTPUT_HC_GVCF}".gz"
tabix -f -p vcf ${OUTPUT_HC_GVCF}".gz"




gatk --java-options "-Xmx48g" HaplotypeCaller  \
   -R ${REF} \
   -L ${INTERVAL} \
   -I ${BAM_PATH} \
   -D ${dbSNP} \
   --dont-use-soft-clipped-bases true \
   -O ${OUTPUT_HC_VCF} 

# 추후 MutationTimeR를 위해서 FORMAT 을 바꿔주기
sed 's/ID=AD,Number=R/ID=AD,Number=2/' ${OUTPUT_HC_VCF} > ${OUTPUT_HC_VCF}".temp"
mv ${OUTPUT_HC_VCF}".temp" ${OUTPUT_HC_VCF}
rm -rf ${OUTPUT_HC_VCF}".temp"

bgzip -c -f ${OUTPUT_HC_VCF} > ${OUTPUT_HC_VCF}".gz"
tabix -f -p vcf ${OUTPUT_HC_VCF}".gz"
