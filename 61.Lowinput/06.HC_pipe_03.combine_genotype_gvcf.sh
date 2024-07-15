#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if ! options=$(getopt -o h --long SCRIPT_COMBINE:,SCRIPT_GENOTYPE:,GENOTYPE_GVCF:,REF:,dbSNP:, -- "$@")
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
        --SCRIPT_COMBINE)
            SCRIPT_COMBINE=$2
        shift 2 ;;
        --SCRIPT_GENOTYPE)
            SCRIPT_GENOTYPE=$2
        shift 2 ;;
        --GENOTYPE_GVCF)
            GENOTYPE_GVCF=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --dbSNP)
            dbSNP=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

# 실행해주기
echo ${SCRIPT_COMBINE}
${SCRIPT_COMBINE}
echo ${SCRIPT_GENOTYPE}
${SCRIPT_GENOTYPE}

gunzip -f -c  ${GENOTYPE_GVCF} > ${GENOTYPE_GVCF%.gz}




################# dbSNP 에 있는건 빼주기 #####################3

GENOTYPE_VCF=${GENOTYPE_GVCF%.gz}

gatk --java-options "-Xmx48g" SelectVariants \
  -R ${REF} \
  -V ${GENOTYPE_VCF} \
  --discordance ${dbSNP} \
  -O ${GENOTYPE_VCF}".filtered"

mv ${GENOTYPE_VCF}".filtered" ${GENOTYPE_VCF}

echo "Subtract dbSNP done"



grep '^#' ${GENOTYPE_VCF} > ${GENOTYPE_VCF}".sorted"
grep -v '^#' ${GENOTYPE_VCF} | sort -k1,1V -k2n >> ${GENOTYPE_VCF}".sorted"
grep -v '^#' ${GENOTYPE_VCF} | grep 'chrX' >> ${GENOTYPE_VCF}".sorted"
grep -v '^#' ${GENOTYPE_VCF} | grep 'chrY' >> ${GENOTYPE_VCF}".sorted"
echo "Sort lexicographically done"

mv ${GENOTYPE_VCF}".sorted" ${GENOTYPE_VCF}
bgzip -c -f ${GENOTYPE_VCF} > ${GENOTYPE_GVCF}
tabix -f -p vcf ${GENOTYPE_GVCF}
echo "bgzip & tabix done"