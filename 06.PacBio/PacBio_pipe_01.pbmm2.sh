#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long BAM_PATH:,REF_MMI:,PRESET:,OUTPUT_BAM:,SAMPLENAME:, -- "$@")
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
        --REF_MMI)
            REF_MMI=$2
        shift 2 ;;
        --PRESET)
            PRESET=$2
        shift 2 ;;
        --OUTPUT_BAM)
            OUTPUT_BAM=$2
        shift 2 ;;
        --SAMPLENAME)
            SAMPLENAME=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

# echo -e "samtools index ${BAM_PATH}"
# samtools index ${BAM_PATH}


echo -e "/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pbmm2 align \
    ${REF_MMI} \
    ${BAM_PATH} \
    ${OUTPUT_BAM} \
    --preset ${PRESET} \
    --sample ${SAMPLENAME}"

# CPU 엄청나게 먹으니 조심해야 할듯
/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pbmm2 align \
    ${REF_MMI} \
    ${BAM_PATH} \
    ${OUTPUT_BAM} \
    --preset ${PRESET} \
    --sample ${SAMPLENAME}

# sort 하면 터진다. samtools로 따로 해주자

mv ${OUTPUT_BAM} ${OUTPUT_BAM%"bam"}"temp.bam"
samtools sort -@ 8 -o ${OUTPUT_BAM} ${OUTPUT_BAM%"bam"}"temp.bam"
samtools index ${OUTPUT_BAM}

rm -rf ${OUTPUT_BAM%"bam"}"temp.bam"