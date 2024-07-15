#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,PRE_BAM_PATH:,MarkDuplicate_PATH:,AddOrReplaceReadGroups_PATH:,BQSR_PATH:,BQSR_RECAL_PATH:,FINAL_BAM_PATH:,REF:,dbsnp:,TMP_PATH:, -- "$@")
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
        --PRE_BAM_PATH)
            PRE_BAM_PATH=$2
        shift 2 ;;
        --MarkDuplicate_PATH)
            MarkDuplicate_PATH=$2
        shift 2 ;;
        --AddOrReplaceReadGroups_PATH)
            AddOrReplaceReadGroups_PATH=$2
        shift 2 ;;
        --BQSR_PATH)
            BQSR_PATH=$2
        shift 2 ;;
        --BQSR_RECAL_PATH)
            BQSR_RECAL_PATH=$2
        shift 2 ;;
        --FINAL_BAM_PATH)
            FINAL_BAM_PATH=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --dbsnp)
            dbsnp=$2
        shift 2 ;;
        --TMP_PATH)
            TMP_PATH=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



# MarkDuplicates & Remove duplicate
java -Xmx60g -jar /opt/Yonsei/Picard/2.26.4/picard.jar MarkDuplicates \
    I=${PRE_BAM_PATH} \
    O=${MarkDuplicate_PATH} \
    M=${MarkDuplicate_PATH%/*}"/"${Sample_ID}'.MarkDuplicate.txt' \
    CREATE_INDEX=true \
    REMOVE_DUPLICATES="true" \
    REMOVE_SEQUENCING_DUPLICATES="true" \
    TMP_DIR=${TMP_PATH}

date
echo "Mark & Remove duplicate done"


# AddOrReplaceRG
java -Xmx60g -jar /opt/Yonsei/Picard/2.26.4/picard.jar AddOrReplaceReadGroups \
    I=${MarkDuplicate_PATH} \
    O=${AddOrReplaceReadGroups_PATH} \
    SORT_ORDER=coordinate \
    RGLB='Meningioma' \
    RGPL='Illumina' \
    RGPU='Illumina' \
    RGSM=${Sample_ID} \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR=${TMP_PATH}

date
echo "AddOrReplaceRG done"



# BaseRecalibrator
    gatk BaseRecalibrator \
    -I ${AddOrReplaceReadGroups_PATH} \
    -R ${REF} \
    --known-sites ${dbsnp} \
    -O ${BQSR_RECAL_PATH} \
    --tmp-dir ${TMP_PATH}

date
echo "1st BaseRecalibrator done"
# ApplyBQSR
    gatk ApplyBQSR \
    -R $REF \
    -bqsr-recal-file ${BQSR_RECAL_PATH} \
    -I ${AddOrReplaceReadGroups_PATH} \
    -O ${BQSR_PATH} \
    --tmp-dir $TMP_PATH
date
echo "1st ApplyBQSR done"


##. Final left Alignment

gatk LeftAlignIndels \
-R $REF \
-I ${BQSR_PATH} \
-O ${FINAL_BAM_PATH}  \
--tmp-dir $TMP_PATH

echo "Left alignment done"
date





