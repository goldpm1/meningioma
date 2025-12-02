#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long OUTPUT_VCF_GZ:,COMMAND:, -- "$@")
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
        --OUTPUT_VCF_GZ)
            OUTPUT_VCF_GZ=$2
        shift 2 ;;
        --COMMAND)
            COMMAND=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


OUTPUT_VCF=${OUTPUT_VCF_GZ%".gz"}

echo -e ${COMMAND}


rm -rf ${OUTPUT_VCF_GZ}  ${OUTPUT_VCF}

${COMMAND}

gunzip ${OUTPUT_VCF_GZ}

# 추후 MutationTimeR를 위해서 FORMAT 을 바꿔주기
sed 's/ID=AD,Number=R/ID=AD,Number=2/' ${OUTPUT_VCF} > ${OUTPUT_VCF}".temp"
mv ${OUTPUT_VCF}".temp" ${OUTPUT_VCF}
rm -rf ${OUTPUT_VCF}".temp" ${OUTPUT_VCF_GZ}

bgzip -c -f ${OUTPUT_VCF} > ${OUTPUT_VCF_GZ}
tabix -f -p vcf ${OUTPUT_VCF_GZ}