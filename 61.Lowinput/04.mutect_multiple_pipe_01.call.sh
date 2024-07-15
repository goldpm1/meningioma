#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long SCRIPT_DIR:,Sample_ID:,BAM_DIR_LIST:,MUTECT_OUTPUT_PATH:,MUTECT_OUTPUT_FMC_PATH:,MUTECT_OUTPUT_FMC_HF_PATH:,MUTECT_OUTPUT_FMC_HF_RMBLACK_PATH:,PON:,REF:,gnomad:,TMP_PATH:,SAMPLE_THRESHOLD:,DP_THRESHOLD:,ALT_THRESHOLD:,REMOVE_MULTIALLELIC:,PASS:,REMOVE_MITOCHONDRIAL_DNA:,BLACKLIST:, -- "$@")
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
        --SCRIPT_DIR)
            SCRIPT_DIR=$2
        shift 2 ;;
        --Sample_ID)
            Sample_ID=$2
        shift 2 ;;
        --BAM_DIR_LIST)
            BAM_DIR_LIST=$2
        shift 2 ;;
        --MUTECT_OUTPUT_PATH)
            MUTECT_OUTPUT_PATH=$2
        shift 2 ;;
        --MUTECT_OUTPUT_FMC_PATH)
            MUTECT_OUTPUT_FMC_PATH=$2
        shift 2 ;;
        --MUTECT_OUTPUT_FMC_HF_PATH)
            MUTECT_OUTPUT_FMC_HF_PATH=$2
        shift 2 ;;
        --MUTECT_OUTPUT_FMC_HF_RMBLACK_PATH)
            MUTECT_OUTPUT_FMC_HF_RMBLACK_PATH=$2
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
        --TMP_PATH)
            TMP_PATH=$2
        shift 2 ;;
        --SAMPLE_THRESHOLD)
            SAMPLE_THRESHOLD=$2
        shift 2 ;;
        --DP_THRESHOLD)
            DP_THRESHOLD=$2
        shift 2 ;;
        --ALT_THRESHOLD)
            ALT_THRESHOLD=$2
        shift 2 ;;
        --REMOVE_MULTIALLELIC)
            REMOVE_MULTIALLELIC=$2
        shift 2 ;;
        --PASS)
            PASS=$2
        shift 2 ;;
        --REMOVE_MITOCHONDRIAL_DNA)
            REMOVE_MITOCHONDRIAL_DNA=$2
        shift 2 ;;
        --BLACKLIST)
            BLACKLIST=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


#1. Multiple (Triple) sample mutect2 call

python3 ${SCRIPT_DIR}"/04.mutect_multiple_pipe_01.call.py" \
    --REF ${REF} \
    --BAM_DIR_LIST ${BAM_DIR_LIST} \
    --normal ${Sample_ID}"_Blood" \
    --panel_of_normals ${PON} \
    --germline_resource ${gnomad} \
    --O ${MUTECT_OUTPUT_PATH} \
    --temp_dir ${TMP_PATH}


# 추후 MutationTimeR를 위해서 FORMAT 을 바꿔주기
sed 's/ID=AD,Number=R/ID=AD,Number=2/' ${MUTECT_OUTPUT_PATH} > ${MUTECT_OUTPUT_PATH}".temp"
mv ${MUTECT_OUTPUT_PATH}".temp" ${MUTECT_OUTPUT_PATH}
rm -rf ${MUTECT_OUTPUT_PATH}".temp"

bgzip -c -f ${MUTECT_OUTPUT_PATH} > ${MUTECT_OUTPUT_PATH}".gz"
tabix -f -p vcf ${MUTECT_OUTPUT_PATH}".gz"


#2. FitlerMutectCall
gatk  FilterMutectCalls -R ${REF} -V ${MUTECT_OUTPUT_PATH} -O ${MUTECT_OUTPUT_FMC_PATH}  \
        --max-events-in-region 1 --min-reads-per-strand 1 --min-median-read-position 8 --min-median-base-quality 20 --min-median-mapping-quality 20

#3. Hard filter
echo -e "python3 mutect_pipe_hardfilter.py --INPUT_VCF "${OUTPUT_FMC_PATH}" --OUTPUT_VCF "${OUTPUT_FMC_HF_PATH}" --SAMPLE_THRESHOLD "${SAMPLE_THRESHOLD}" --DP_THRESHOLD "${DP_THRESHOLD}" --ALT_THRESHOLD "${ALT_THRESHOLD}"  --REMOVE_MULTIALLELIC "${REMOVE_MULTIALLELIC}" --PASS "${PASS}" --REMOVE_MITOCHONDRIAL_DNA "${REMOVE_MITOCHONDRIAL_DNA}
python3 /data/project/Meningioma/script/03.Variant_calling\&Annotation/mutect_pipe_hardfilter.py \
        --INPUT_VCF ${MUTECT_OUTPUT_FMC_PATH} \
        --OUTPUT_VCF ${MUTECT_OUTPUT_FMC_HF_PATH} \
        --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD}  --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA}

#4. Repeat region 지우기
bedtools intersect -header -v -a ${MUTECT_OUTPUT_FMC_HF_PATH} -b ${BLACKLIST} > ${MUTECT_OUTPUT_FMC_HF_RMBLACK_PATH}
