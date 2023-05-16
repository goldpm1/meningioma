#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long REF:,DOC_PATH:,FINAL_BAM_PATH:,INTERVAL:, -- "$@")
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
        --DOC_PATH)
            DOC_PATH=$2
        shift 2 ;;
        --FINAL_BAM_PATH)
            FINAL_BAM_PATH=$2
        shift 2 ;;
        --INTERVAL)
            INTERVAL=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo -e "java -Xmx12g -jar /opt/Yonsei/GATK/3.8-1/GenomeAnalysisTK.jar  -T DepthOfCoverage  -R "${REF}" -o "${DOC_PATH}" -I "${FINAL_BAM_PATH}" -L "${INTERVAL}


# DepthOfCoverage #
# java -Xmx12g -jar /opt/Yonsei/GATK/3.8-1/GenomeAnalysisTK.jar  -T DepthOfCoverage \
#     -R ${REF} \
#     -o ${DOC_PATH} \
#     -I ${FINAL_BAM_PATH} \
#     -L ${INTERVAL}
    
gatk DepthOfCoverage \
    -R ${REF} \
    -O ${DOC_PATH} \
    -I ${FINAL_BAM_PATH} \
    -L ${INTERVAL}
    
#/data/project/Alzheimer/UCSC.hg38.WholeGenome.bed


date
echo "DepthofCoverage done"
