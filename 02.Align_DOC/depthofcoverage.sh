#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long REF:,DOC_PATH:,BAM_PATH:,INTERVAL:, -- "$@")
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
        --BAM_PATH)
            BAM_PATH=$2
        shift 2 ;;
        --INTERVAL)
            INTERVAL=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo -e "java -Xmx12g -jar /opt/Yonsei/GATK/3.8-1/GenomeAnalysisTK.jar  -T DepthOfCoverage  -R "${REF}" -o "${DOC_PATH}" -I "${BAM_PATH}" -L "${INTERVAL}
java -Xmx12g -jar /opt/Yonsei/GATK/3.8-1/GenomeAnalysisTK.jar  -T DepthOfCoverage \
    -R ${REF} \
    -o ${DOC_PATH} \
    -I ${BAM_PATH} \
    -L ${INTERVremote

# echo -e "java -Xmx12g -jar /opt/Yonsei/GATK/4.2.3.0/gatk-package-4.2.3.0-local.jar DepthOfCoverage  -R ${REF}  -O ${DOC_PATH}  -I ${BAM_PATH}  -L ${INTERVAL}"
# java -Xmx12g -jar /opt/Yonsei/GATK/4.2.3.0/gatk-package-4.2.3.0-local.jar DepthOfCoverage \
#     -R ${REF} \
#     -O ${DOC_PATH} \
#     -I ${BAM_PATH} \
#     -L ${INTERVAL} \
#     --disable-tool-default-read-filters

#echo -e "gatk DepthOfCoverage  -R "${REF}" -O "${DOC_PATH}" -I "${BAM_PATH}" -L "${INTERVAL}
# gatk DepthOfCoverage \
#     -R ${REF} \
#     -O ${DOC_PATH} \
#     -I ${BAM_PATH} \
#     -L ${INTERVAL} 

    
#/data/project/Alzheimer/UCSC.hg38.WholeGenome.bed


date
echo "DepthofCoverage done"
