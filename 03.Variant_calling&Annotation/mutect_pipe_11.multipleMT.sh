#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,CASE_BAM_PATH1:,CASE_BAM_PATH2:,CONTROL_BAM_PATH:,OUTPUT_VCF_GZ:,OUTPUT_FMC_PATH:,OUTPUT_FMC_HF_PATH:,OUTPUT_FMC_HF_RMBLACK_PATH:,PON:,REF:,gnomad:,INTERVAL:,TMP_PATH:,SAMPLE_THRESHOLD:,DP_THRESHOLD:,ALT_THRESHOLD:,REMOVE_MULTIALLELIC:,PASS:,REMOVE_MITOCHONDRIAL_DNA:,BLACKLIST:, -- "$@")
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
        --CASE_BAM_PATH1)
            CASE_BAM_PATH1=$2
        shift 2 ;;
        --CASE_BAM_PATH2)
            CASE_BAM_PATH2=$2
        shift 2 ;;
        --CONTROL_BAM_PATH)
            CONTROL_BAM_PATH=$2
        shift 2 ;;
        --OUTPUT_VCF_GZ)
            OUTPUT_VCF_GZ=$2
        shift 2 ;;
        --OUTPUT_FMC_PATH)
            OUTPUT_FMC_PATH=$2
        shift 2 ;;
        --OUTPUT_FMC_HF_PATH)
            OUTPUT_FMC_HF_PATH=$2
        shift 2 ;;
        --OUTPUT_FMC_HF_RMBLACK_PATH)
            OUTPUT_FMC_HF_RMBLACK_PATH=$2
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
        --INTERVAL)
            INTERVAL=$2
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


#1. Multiple sample mutect2 call
gatk --java-options "-Xmx48g" Mutect2 \
-R $REF \
-I ${CASE_BAM_PATH1} \
-I ${CASE_BAM_PATH2} \
-I ${CONTROL_BAM_PATH} \
-normal ${Sample_ID}"_Blood" \
--panel-of-normals ${PON} \
--germline-resource ${gnomad} \
-O ${OUTPUT_VCF_GZ} \
--tmp-dir ${TMP_PATH}


gunzip ${OUTPUT_VCF_GZ}
OUTPUT_VCF=${OUTPUT_VCF_GZ%".gz"}

# 추후 MutationTimeR를 위해서 FORMAT 을 바꿔주기
sed 's/ID=AD,Number=R/ID=AD,Number=2/' ${OUTPUT_VCF} > ${OUTPUT_VCF}".temp"
mv ${OUTPUT_VCF}".temp" ${OUTPUT_VCF}
rm -rf ${OUTPUT_VCF}".temp"

bgzip -c -f ${OUTPUT_VCF} > ${OUTPUT_VCF_GZ}
tabix -p vcf ${OUTPUT_VCF_GZ}



#2. FitlerMutectCall
gatk  FilterMutectCalls -R ${REF} -V ${OUTPUT_VCF_GZ} -O ${OUTPUT_FMC_PATH} --max-events-in-region 1 --min-reads-per-strand 1 --min-median-read-position 8 --min-median-base-quality 20 --min-median-mapping-quality 20

# 추후 MutationTimeR를 위해서 FORMAT 을 바꿔주기
sed 's/ID=AD,Number=R/ID=AD,Number=2/' ${OUTPUT_FMC_PATH} > ${OUTPUT_FMC_PATH}".temp"
mv ${OUTPUT_FMC_PATH}".temp" ${OUTPUT_FMC_PATH}
rm -rf ${OUTPUT_FMC_PATH}".temp"


#3. Hard filter
echo -e "python3 mutect_pipe_hardfilter.py --INPUT_VCF "${OUTPUT_FMC_PATH}" --OUTPUT_VCF "${OUTPUT_FMC_HF_PATH}" --SAMPLE_THRESHOLD "${SAMPLE_THRESHOLD}" --DP_THRESHOLD "${DP_THRESHOLD}" --ALT_THRESHOLD "${ALT_THRESHOLD}"  --REMOVE_MULTIALLELIC "${REMOVE_MULTIALLELIC}" --PASS "${PASS}" --REMOVE_MITOCHONDRIAL_DNA "${REMOVE_MITOCHONDRIAL_DNA}
python3 mutect_pipe_hardfilter.py --INPUT_VCF ${OUTPUT_FMC_PATH} --OUTPUT_VCF ${OUTPUT_FMC_HF_PATH} --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD}  --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA}

#4. Repeat region 지우기
bedtools intersect -header -v -a ${OUTPUT_FMC_HF_PATH} -b ${BLACKLIST} > ${OUTPUT_FMC_HF_RMBLACK_PATH}











# bgzip -c -f ${OUTPUT_FMC_PATH%".gz"} > ${OUTPUT_FMC_PATH}
# tabix -p vcf ${OUTPUT_FMC_PATH
