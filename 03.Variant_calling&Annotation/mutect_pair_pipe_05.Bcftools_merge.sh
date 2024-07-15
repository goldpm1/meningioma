#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long BCFTOOLS_MERGE_TXT:,BCFTOOLS_MERGE_VCF_GZ:,BCFTOOLS_MERGE_VCF:, -- "$@")
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
        --BCFTOOLS_MERGE_TXT)
            BCFTOOLS_MERGE_TXT=$2
        shift 2 ;;
        --BCFTOOLS_MERGE_VCF_GZ)
            BCFTOOLS_MERGE_VCF_GZ=$2
        shift 2 ;;
        --BCFTOOLS_MERGE_VCF)
            BCFTOOLS_MERGE_VCF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo -e "bcftools merge --force-samples --merge none -l ${BCFTOOLS_MERGE_TXT} -O v -o ${BCFTOOLS_MERGE_VCF}"

bcftools merge --force-samples --merge none -l ${BCFTOOLS_MERGE_TXT} -O v -o ${BCFTOOLS_MERGE_VCF}

# sort lexicographically
grep "^#" ${BCFTOOLS_MERGE_VCF} > ${BCFTOOLS_MERGE_VCF}".sorted" && grep -v "^#" ${BCFTOOLS_MERGE_VCF} | sort -k1,1V -k2n >> ${BCFTOOLS_MERGE_VCF}".sorted"
mv ${BCFTOOLS_MERGE_VCF}".sorted" ${BCFTOOLS_MERGE_VCF}
rm -rf ${BCFTOOLS_MERGE_VCF}".sorted"

# bgzip & tabix
bgzip -c -f ${BCFTOOLS_MERGE_VCF} > ${BCFTOOLS_MERGE_VCF_GZ}
tabix -p vcf ${BCFTOOLS_MERGE_VCF_GZ}