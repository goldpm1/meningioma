#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,CASE_BAM_PATH:,CONTROL_BAM_PATH:,OUTPUT_VCF_GZ:,PON:,REF:,gnomad:,INTERVAL:,TMP_PATH:, -- "$@")
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
        --Sample_ID)
            Sample_ID=$2
        shift 2 ;;
        --CASE_BAM_PATH)
            CASE_BAM_PATH=$2
        shift 2 ;;
        --CONTROL_BAM_PATH)
            CONTROL_BAM_PATH=$2
        shift 2 ;;
        --OUTPUT_VCF_GZ)
            OUTPUT_VCF_GZ=$2
        shift 2 ;;
        --PON)
            PON=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --gnomad)
            gnomad=$2
        shift 2 ;;
        --INTERVAL)
            INTERVAL=$2
        shift 2 ;;
        --TMP_PATH)
            TMP_PATH=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo $CONTROL_BAM_PATH"/Blood_WGS.bam"


OUTPUT_VCF=${OUTPUT_VCF_GZ%".gz"}

rm -rf ${OUTPUT_VCF_GZ}  ${OUTPUT_VCF}

#1.Paired sample
gatk --java-options "-Xmx48g" Mutect2 \
-R $REF \
-I ${CASE_BAM_PATH} \
-I ${CONTROL_BAM_PATH} \
-normal ${Sample_ID}"_Blood" \
--panel-of-normals $PON \
--germline-resource $gnomad \
-O ${OUTPUT_VCF_GZ} \
--tmp-dir $TMP_PATH \
--intervals $INTERVAL  \
--dont-use-soft-clipped-bases true 


gunzip ${OUTPUT_VCF_GZ}

# 추후 MutationTimeR를 위해서 FORMAT 을 바꿔주기
sed 's/ID=AD,Number=R/ID=AD,Number=2/' ${OUTPUT_VCF} > ${OUTPUT_VCF}".temp"
mv ${OUTPUT_VCF}".temp" ${OUTPUT_VCF}
rm -rf ${OUTPUT_VCF}".temp" ${OUTPUT_VCF_GZ}

bgzip -c -f ${OUTPUT_VCF} > ${OUTPUT_VCF_GZ}
tabix -f -p vcf ${OUTPUT_VCF_GZ}