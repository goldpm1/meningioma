#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long REF:,REF_MMI:,PRESET:, -- "$@")
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
        --REF_MMI)
            REF_MMI=$2
        shift 2 ;;
        --PRESET)
            PRESET=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

# echo -e "samtools index ${BAM_PATH}"
# samtools index ${BAM_PATH}


echo -e "/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pbmm2 index \
    ${REF} \
    ${REF_MMI} \
    --preset ${PRESET}"

/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pbmm2 index \
    ${REF} \
    ${REF_MMI} \
    --preset ${PRESET}
