#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long BAM_PATH:,OUTPUT_SIG:,REF:,OUTPUT_VCF:, -- "$@")
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
        --OUTPUT_SIG)
            OUTPUT_SIG=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --OUTPUT_VCF)
            OUTPUT_VCF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo -e "/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pbsv discover \
    --hifi \
    --region "chr22" \
    ${BAM_PATH} \
    ${OUTPUT_SIG}"

/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pbsv discover \
    --hifi \
    --region "chr22" \
    ${BAM_PATH} \
    ${OUTPUT_SIG}



echo -e "/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pbsv call \
    --hifi \
    --region "chr22" \
    ${REF} \
    ${OUTPUT_SIG} \
    ${OUTPUT_VCF}"

/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pbsv call \
    --hifi \
    --region "chr22" \
    ${REF} \
    ${OUTPUT_SIG} \
    ${OUTPUT_VCF}

