#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if ! options=$(getopt -o h --long BAM_PATH:,INTERVAL:,REF:,dbSNP:,OUTPUT_HC:, -- "$@")
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
        --OUTPUT_HC)
            OUTPUT_HC=$2
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
   -O ${OUTPUT_HC} 

# 추후 MutationTimeR를 위해서 FORMAT 을 바꿔주기
sed 's/ID=AD,Number=R/ID=AD,Number=2/' ${OUTPUT_HC} > ${OUTPUT_HC}".temp"
mv ${OUTPUT_HC}".temp" ${OUTPUT_HC}
rm -rf ${OUTPUT_HC}".temp"

bgzip -c -f ${OUTPUT_HC} > ${OUTPUT_HC}".gz"
tabix -f -p vcf ${OUTPUT_HC}".gz"