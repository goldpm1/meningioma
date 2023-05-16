#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,CASE_BAM_PATH:,CONTROL_BAM_PATH:,OUTPUT_VCF_GZ:,OUTPUT_FMC_PATH:,OUTPUT_FMC_HF_PATH:,OUTPUT_FMC_HF_RMBLACK_PATH:,PON:,REF:,gnomad:,INTERVAL:,TMP_PATH:,SAMPLE_THRESHOLD:,DP_THRESHOLD:,ALT_THRESHOLD:,REMOVE_MULTIALLELIC:,PASS:,REMOVE_MITOCHONDRIAL_DNA:,BLACKLIST:, -- "$@")
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
        --CASE_BAM_PATH)
            CASE_BAM_PATH=$2
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



echo -e "REF : "${REF}" OUTPUT_VCF_GZ :"${OUTPUT_VCF_GZ}" OUTPUT_FMC_PATH : "${OUTPUT_FMC_PATH}

#2. FilterMutectCall
gatk FilterMutectCalls -R ${REF} -V ${OUTPUT_VCF_GZ} -O ${OUTPUT_FMC_PATH} --max-events-in-region 1 --min-median-read-position 8 --min-median-base-quality 20  --min-reads-per-strand 1 --min-median-mapping-quality 20

# 추후 MutationTimeR를 위해서 FORMAT 을 바꿔주기
sed 's/ID=AD,Number=R/ID=AD,Number=2/' ${OUTPUT_FMC_PATH} > ${OUTPUT_FMC_PATH}".temp"
mv ${OUTPUT_FMC_PATH}".temp" ${OUTPUT_FMC_PATH}
rm -rf ${OUTPUT_FMC_PATH}".temp"


 # min-reads-per-strand forward, reverse strand가 적어도 1개씩은 있어야 뽑느다

#3. Hard filter
echo -e "python3 mutect_pipe_hardfilter.py --INPUT_VCF "${OUTPUT_FMC_PATH}" --OUTPUT_VCF "${OUTPUT_FMC_HF_PATH}" --SAMPLE_THRESHOLD "${SAMPLE_THRESHOLD}" --DP_THRESHOLD "${DP_THRESHOLD}" --ALT_THRESHOLD "${ALT_THRESHOLD}"  --REMOVE_MULTIALLELIC "${REMOVE_MULTIALLELIC}" --PASS "${PASS}" --REMOVE_MITOCHONDRIAL_DNA "${REMOVE_MITOCHONDRIAL_DNA}
python3 mutect_pipe_hardfilter.py --INPUT_VCF ${OUTPUT_FMC_PATH} --OUTPUT_VCF ${OUTPUT_FMC_HF_PATH} --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD}  --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA}

#4. Repeat region 지우기
bedtools intersect -header -v -a ${OUTPUT_FMC_HF_PATH} -b ${BLACKLIST} > ${OUTPUT_FMC_HF_RMBLACK_PATH}











# bgzip -c -f ${OUTPUT_FMC_PATH%".gz"} > ${OUTPUT_FMC_PATH}
# tabix -p vcf ${OUTPUT_FMC_PATH
