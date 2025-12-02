#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if ! options=$(getopt -o h --long INPUT_GVCF:,OUTPUT_GVCF:, -- "$@")
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
        --INPUT_GVCF)
            INPUT_GVCF=$2
        shift 2 ;;
        --OUTPUT_GVCF)
            OUTPUT_GVCF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

# <NON_REF> 단독은 없애주기
grep '#' ${INPUT_GVCF} > ${OUTPUT_GVCF%g.vcf}"pre1.g.header.vcf"
grep -v '#' ${INPUT_GVCF}  | awk '{ if ($5 != "<NON_REF>") { print } }'  > ${OUTPUT_GVCF%g.vcf}"pre1.g.body.vcf"
cat ${OUTPUT_GVCF%g.vcf}"pre1.g.header.vcf" ${OUTPUT_GVCF%g.vcf}"pre1.g.body.vcf" > ${OUTPUT_GVCF%g.vcf}"PASS.pre1.g.vcf"


# biallelic을 monoallelic로 변환 후 NON_REF를 없애주기
bcftools norm -m - ${OUTPUT_GVCF%g.vcf}"PASS.pre1.g.vcf"  -O v -o ${OUTPUT_GVCF%g.vcf}"PASS.pre2.g.vcf"
grep '#'  ${OUTPUT_GVCF%g.vcf}"PASS.pre1.g.vcf" > ${OUTPUT_GVCF%g.vcf}"PASS.pre3.g.vcf"
grep -v '#' ${OUTPUT_GVCF%g.vcf}"PASS.pre2.g.vcf"  | awk '{ if ($5 != "<NON_REF>") { print } }'  >> ${OUTPUT_GVCF%g.vcf}"PASS.pre3.g.vcf"


# 다시 multiallelic으로 합쳐주기
bcftools norm -m +  ${OUTPUT_GVCF%g.vcf}"PASS.pre3.g.vcf" >  ${OUTPUT_GVCF%g.vcf}"PASS.pre4.g.vcf"
bgzip -c  ${OUTPUT_GVCF%g.vcf}"PASS.pre4.g.vcf" >  ${OUTPUT_GVCF%g.vcf}"PASS.pre4.g.vcf.gz"
tabix -p vcf ${OUTPUT_GVCF%g.vcf}"PASS.pre4.g.vcf.gz"

# multiallelic 제거하기
bcftools view --max-alleles 2 ${OUTPUT_GVCF%g.vcf}"PASS.pre4.g.vcf.gz" > ${OUTPUT_GVCF}


bgzip -c -f ${OUTPUT_GVCF} > ${OUTPUT_GVCF}".gz"
tabix -p vcf ${OUTPUT_GVCF}".gz"
#tabix -p vcf ${OUTPUT_GVCF}

rm -rf ${OUTPUT_GVCF%g.vcf}"pre1.g.header.vcf" ${OUTPUT_GVCF%g.vcf}"pre1.g.body.vcf" ${OUTPUT_GVCF%g.vcf}"PASS.pre1.g.vcf" ${OUTPUT_GVCF%g.vcf}"PASS.pre2.g.vcf" ${OUTPUT_GVCF%g.vcf}"PASS.pre3.g.vcf" ${OUTPUT_GVCF%g.vcf}"PASS.pre4.g.vcf" ${OUTPUT_GVCF%g.vcf}"PASS.pre4.g.vcf.gz"  ${OUTPUT_GVCF%g.vcf}"PASS.pre4.g.vcf.gz.tbi"