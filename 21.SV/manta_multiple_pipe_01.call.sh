#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long CONTROL_BAM_PATH:,CASE_BAM_PATH1:,CASE_BAM_PATH2:,REF:,CALLREGIONS:,OUTPUT_DIR:,OUTPUT_PATH_SOMATIC:,OUTPUT_PATH_SOMATIC_GZ:,OUTPUT_PATH_DIPLOID:,OUTPUT_PATH_DIPLOID_GZ:, -- "$@")
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
        --CONTROL_BAM_PATH)
            CONTROL_BAM_PATH=$2
        shift 2 ;;
        --CASE_BAM_PATH1)
            CASE_BAM_PATH1=$2
        shift 2 ;;
        --CASE_BAM_PATH2)
            CASE_BAM_PATH2=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --CALLREGIONS)
            CALLREGIONS=$2
        shift 2 ;;
        --OUTPUT_DIR)
            OUTPUT_DIR=$2
        shift 2 ;;
        --OUTPUT_PATH_SOMATIC)
            OUTPUT_PATH_SOMATIC=$2
        shift 2 ;;
        --OUTPUT_PATH_SOMATIC_GZ)
            OUTPUT_PATH_SOMATIC_GZ=$2
        shift 2 ;;
        --OUTPUT_PATH_DIPLOID)
            OUTPUT_PATH_DIPLOID=$2
        shift 2 ;;
        --OUTPUT_PATH_DIPLOID_GZ)
            OUTPUT_PATH_DIPLOID_GZ=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


# #1. Configuration

echo -e "/opt/Yonsei/python/2.7.17/bin/python2 /opt/Yonsei/Manta/1.6.0/bin/configManta.py \
--referenceFasta $REF \
--normalBam ${CONTROL_BAM_PATH} \
--tumorBam ${CASE_BAM_PATH1}  \
--runDir ${OUTPUT_DIR} \
--generateEvidenceBam \
--exome \
--callRegions ${CALLREGIONS} "

/opt/Yonsei/python/2.7.17/bin/python2 /opt/Yonsei/Manta/1.6.0/bin/configManta.py \
--referenceFasta $REF \
--normalBam ${CONTROL_BAM_PATH} \
--tumorBam ${CASE_BAM_PATH1}  \
--runDir ${OUTPUT_DIR} \
--generateEvidenceBam \
--exome \
--callRegions ${CALLREGIONS}

date
echo "Configuration done"


# #2. Execution
/opt/Yonsei/python/2.7.17/bin/python2 ${OUTPUT_DIR}"/runWorkflow.py" -j 12
date
echo "Execution done"

#tumor only는 somaticSV,  공통으로 다 나온거면 diploidSV
zcat ${OUTPUT_DIR}"/results/variants/somaticSV.vcf.gz" | grep '#' > ${OUTPUT_PATH_SOMATIC}   # 일단 header만 복사
zcat ${OUTPUT_DIR}"/results/variants/somaticSV.vcf.gz" | grep -v 'IMPRECISE' | awk -F "\t" '{if($7=="PASS"){print $0}}' >> ${OUTPUT_PATH_SOMATIC}
bgzip -c -f  ${OUTPUT_PATH_SOMATIC} > ${OUTPUT_PATH_SOMATIC_GZ}
tabix -f -p vcf ${OUTPUT_PATH_SOMATIC_GZ}


zcat ${OUTPUT_DIR}"/results/variants/diploidSV.vcf.gz" | grep '#' > ${OUTPUT_PATH_DIPLOID}   # 일단 header만 복사
zcat ${OUTPUT_DIR}"/results/variants/diploidSV.vcf.gz" | grep -v 'IMPRECISE' | awk -F "\t" '{if($7=="PASS"){print $0}}' >> ${OUTPUT_PATH_DIPLOID}
bgzip -c -f  ${OUTPUT_PATH_DIPLOID} > ${OUTPUT_PATH_DIPLOID_GZ}
tabix -f -p vcf ${OUTPUT_PATH_DIPLOID_GZ}