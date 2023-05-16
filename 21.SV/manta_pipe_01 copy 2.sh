#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long CONTROL_BAM_PATH:,CASE_BAM_PATH:,REF:,CALLREGIONS:,OUTPUT_PATH_01:,OUTPUT_PATH_02:, -- "$@")
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
        --CASE_BAM_PATH)
            CASE_BAM_PATH=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --CALLREGIONS)
            CALLREGIONS=$2
        shift 2 ;;
        --OUTPUT_PATH_01)
            OUTPUT_PATH_01=$2
        shift 2 ;;
        --OUTPUT_PATH_02)
            OUTPUT_PATH_02=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


# #1. Configuration

echo -e "/opt/Yonsei/python/2.7.17/bin/python2 /opt/Yonsei/Manta/1.6.0/bin/configManta.py \
--referenceFasta $REF \
--tumorBam ${CASE_BAM_PATH} \
--normalBam ${CONTROL_BAM_PATH}  \
--runDir ${OUTPUT_PATH_01} \
--generateEvidenceBam \
--exome \
--callRegions ${CALLREGIONS} "

/opt/Yonsei/python/2.7.17/bin/python2 /opt/Yonsei/Manta/1.6.0/bin/configManta.py \
--referenceFasta $REF \
--tumorBam ${CASE_BAM_PATH} \
--normalBam ${CONTROL_BAM_PATH}  \
--runDir ${OUTPUT_PATH_01} \
--generateEvidenceBam \
--exome \
--callRegions ${CALLREGIONS}

date
echo "Configuration done"


# #2. Execution
/opt/Yonsei/python/2.7.17/bin/python2 ${OUTPUT_PATH_01}"/runWorkflow.py" -j 12
date
echo "Execution done"

#tumor only는 somaticSV,  둘 다 나온거면 diploidSV
zcat ${OUTPUT_PATH_01}"/results/variants/diploidSV.vcf.gz" | grep '#' > ${OUTPUT_PATH_02}   # 일단 header만 복사
zcat ${OUTPUT_PATH_01}"/results/variants/diploidSV.vcf.gz" | grep -v 'IMPRECISE' | awk -F "\t" '{if($7=="PASS"){print $0}}' >> ${OUTPUT_PATH_02}

bgzip -c -f  ${OUTPUT_PATH_02} > ${OUTPUT_PATH_02}".gz"
tabix -p vcf ${OUTPUT_PATH_02}".gz"