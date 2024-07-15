#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if ! options=$(getopt -o h --long REF:,dbSNP:,BLACKLIST:,OUTPUT_HC:,OUTPUT_HC_VF:, -- "$@")
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
        --REF)
            REF=$2
        shift 2 ;;
        --dbSNP)
            dbSNP=$2
        shift 2 ;;
        --OUTPUT_HC)
            OUTPUT_HC=$2
        shift 2 ;;
        --BLACKLIST)
            BLACKLIST=$2
        shift 2 ;;
        --OUTPUT_HC_VF)
            OUTPUT_HC_VF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

rm -rf ${OUTPUT_HC_VF}

#확장자
END=${OUTPUT_HC##*.}
echo ${END}

gatk --java-options "-Xmx48g" VariantFiltration  \
   -R ${REF} \
   -V ${OUTPUT_HC} \
   -O ${OUTPUT_HC_VF%.*}".VF."${END} \
   --filter-name "my_filter1" \
   --filter-expression "DP<=24"
echo -e "VariantFiltration done\n"

if [ "${dbSNP}" == "False" ]; then
    gatk --java-options "-Xmx48g" SelectVariants \
    -R ${REF} \
    -V ${OUTPUT_HC_VF%.*}".VF."${END} \
    -O ${OUTPUT_HC_VF%.*}".VF.SV."${END}
else
    gatk --java-options "-Xmx48g" SelectVariants \
    -R ${REF} \
    -V ${OUTPUT_HC_VF%.*}".VF."${END} \
    --discordance ${dbSNP} \
    -O ${OUTPUT_HC_VF%.*}".VF.SV."${END}
fi
echo -e "SelectVariants done\n"

grep '#' ${OUTPUT_HC_VF%.*}".VF.SV."${END} > ${OUTPUT_HC_VF} 
grep 'PASS' ${OUTPUT_HC_VF%.*}".VF.SV."${END} |  grep -v 'chrM' >> ${OUTPUT_HC_VF} 
rm -rf ${OUTPUT_HC_VF%.*}".VF."${END}  ${OUTPUT_HC_VF%.*}".VF."${END}".idx" ${OUTPUT_HC_VF%.*}".VF.SV."${END} ${OUTPUT_HC_VF%.*}".VF.SV."${END}".idx"

echo -e "Select PASS done\n"


bedtools intersect -header -v -a ${OUTPUT_HC_VF} -b ${BLACKLIST} > ${OUTPUT_HC_VF}".rmblack.vcf"
mv ${OUTPUT_HC_VF}".rmblack.vcf" ${OUTPUT_HC_VF}
rm -rf ${OUTPUT_HC_VF}".rmblack.vcf"  ${OUTPUT_HC_VF}"temp2.vcf" 

echo -e "Remove blacklist done\n"

bgzip -c -f ${OUTPUT_HC_VF} > ${OUTPUT_HC_VF}".gz"
tabix -f -p vcf ${OUTPUT_HC_VF}".gz"


# bash "/data/project/OPLL/script/43.Variant_calling_annotation/mutect_pipe_20.vep.sh" \
#     --REF ${REF} \
#     --INPUT_VCF  ${OUTPUT_HC_VF} \
#     --OUTPUT_VCF ${OUTPUT_HC_VF_VEP}