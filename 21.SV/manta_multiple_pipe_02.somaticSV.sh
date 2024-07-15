#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long INPUT_PATH:,OUTPUT_PATH:,OUTPUT_INV_PATH:,REF:,PANDAS_DIR:,ID:,CONTROL_SAMPLE:, -- "$@")
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
        --INPUT_PATH)
            INPUT_PATH=$2
        shift 2 ;;
        --OUTPUT_PATH)
            OUTPUT_PATH=$2
        shift 2 ;;
        --OUTPUT_INV_PATH)
            OUTPUT_INV_PATH=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --PANDAS_DIR)
            PANDAS_DIR=$2
        shift 2 ;;
        --ID)
            ID=$2
        shift 2 ;;
        --CONTROL_SAMPLE)
            CONTROL_SAMPLE=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



# [Pandas parsing]
python3 "manta_multiple_somatic_vcftobed.py" \
    --INPUT_PATH ${INPUT_PATH} \
    --OUTPUT_PATH ${OUTPUT_PATH} \
    --PANDAS_DIR ${PANDAS_DIR} \
    --ID ${ID} \
    --CONTROL_SAMPLE ${CONTROL_SAMPLE}
echo "Pandas parsing done"

# [convertInversion]
python /data/project/DC_WGS/script/convertInversion.py \
/opt/Yonsei/samtools/1.7/samtools \
${REF} \
${OUTPUT_PATH} \
> ${OUTPUT_INV_PATH}
echo "convertInversion done"

# [Pandas parsing : INV]
python3 manta_vcftobed_inv.py \
    --INPUT_PATH ${OUTPUT_INV_PATH} \
    --PANDAS_DIR  ${PANDAS_DIR} \
    --ID ${ID} \
    --CONTROL_SAMPLE ${CONTROL_SAMPLE}

rm -rf ${OUTPUT_PATH}
echo "Pandas_inv parsing done"


# [Pandas parsing한 것 2개 합치고 옮겨주기]
# cat $OutPath"/3.pandas/"$ID"/"$ID".deldup.Manta.bed" $OutPath"/3.pandas/"$ID"/"$ID".inv.Manta.bed" | sort -V -k1,1 -k2,2 \
# > $OutPath"/3.pandas/"$ID"/"$ID".Manta.bed"

# [Bedpe 파일 만들기 : 양쪽에 길이의 10%정도만큼 늘려준다 ]
#python3 bedtobedpe.py --InputBed $OutPath"/3.pandas/"$ID"/"$ID".Manta.bed"  --ExtensionFrac 0.1
