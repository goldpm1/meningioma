#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long Sample_ID:,TISSUE_BROAD:,OUTPUT_VCF_GZ:,OUTPUT_FMC_PATH:,OUTPUT_FMC_HF_PATH:,OUTPUT_FMC_HF_RMBLACK_PATH:,PON:,REF:,gnomad:,INTERVAL:,TMP_PATH:,SAMPLE_THRESHOLD:,DP_THRESHOLD:,ALT_THRESHOLD:,TLOD_THRESHOLD:,MIN_BQ:,REMOVE_MULTIALLELIC:,PASS:,REMOVE_MITOCHONDRIAL_DNA:,BLACKLIST:, -- "$@")
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
        --TISSUE_BROAD)
            TISSUE_BROAD=$2
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
        --TLOD_THRESHOLD)
            TLOD_THRESHOLD=$2
        shift 2 ;;
        --MIN_BQ)
            MIN_BQ=$2
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


echo "Start time: $(date)\n"
echo -e "REF : "${REF}"\nOUTPUT_VCF_GZ :"${OUTPUT_VCF_GZ}"\nOUTPUT_FMC_PATH : "${OUTPUT_FMC_PATH}"\n\n"

if [ ! -d ${OUTPUT_FMC_PATH%/*} ] ; then
    mkdir ${OUTPUT_FMC_PATH%/*}
fi


#. FilterMutectCall
gatk FilterMutectCalls -R ${REF} -V ${OUTPUT_VCF_GZ} -O ${OUTPUT_FMC_PATH} --max-events-in-region 1 --min-median-read-position 8 --min-median-base-quality ${MIN_BQ}  --min-reads-per-strand 1 --min-median-mapping-quality 20

# VariantFiltration by TLOD_THRESHOLD
# "TLOD_THRESHOLD" 값이 0이 아니고, 비어있지 않은 경우에만 아래의 블록을 실행한다는 뜻
if [[ "${TLOD_THRESHOLD}" != "0" && -n "${TLOD_THRESHOLD}" ]]; then
    echo "TLOD filter ( < ${TLOD_THRESHOLD}) start: $(date)"
    bcftools norm -m -any -Oz -o ${OUTPUT_FMC_PATH%vcf}"split.vcf.gz" ${OUTPUT_FMC_PATH}
    gatk IndexFeatureFile -I ${OUTPUT_FMC_PATH%vcf}"split.vcf.gz"
    gatk VariantFiltration \
        -R ${REF} \
        -V ${OUTPUT_FMC_PATH%vcf}"split.vcf.gz" \
        -O ${OUTPUT_FMC_PATH%vcf}"TLODfilter.vcf" \
        --filter-expression "INFO.TLOD < ${TLOD_THRESHOLD}" \
        --filter-name "Low_TLOD"

    rm -rf ${OUTPUT_FMC_PATH%vcf}"split.vcf.gz" ${OUTPUT_FMC_PATH%vcf}"split.vcf.gz.tbi" # 만들어놨던 중간 파일 삭제
    cp ${OUTPUT_FMC_PATH%vcf}"TLODfilter.vcf" ${OUTPUT_FMC_PATH}
fi

# 추후 MutationTimeR를 위해서 FORMAT 을 바꿔주기
echo -e "sed -i '3i ##FILTER=<ID=RESCUE,Description="Rescued by Pysam">' ${OUTPUT_FMC_PATH}"
sed -i '3i ##FILTER=<ID=RESCUE,Description="Rescued by Pysam">' ${OUTPUT_FMC_PATH}
sed 's/ID=AD,Number=R/ID=AD,Number=2/' ${OUTPUT_FMC_PATH} > ${OUTPUT_FMC_PATH}".temp"
mv ${OUTPUT_FMC_PATH}".temp" ${OUTPUT_FMC_PATH}
rm -rf ${OUTPUT_FMC_PATH}".temp"


# chr22:29658232 위치의 FILTER를 PASS로 변경
awk 'BEGIN{FS=OFS="\t"} $1=="chr22" && $2=="29658232" {$7="PASS"} {print}' ${OUTPUT_FMC_PATH} > ${OUTPUT_FMC_PATH}".temp"
mv ${OUTPUT_FMC_PATH}".temp" ${OUTPUT_FMC_PATH}



#3. Hard filter ("PASS"만 선택하기)
echo "Hard filter start: $(date)"
echo -e "python3 mutect_pipe_hardfilter.py --INPUT_VCF "${OUTPUT_FMC_PATH}" --OUTPUT_VCF "${OUTPUT_FMC_HF_PATH}" --SAMPLE_THRESHOLD "${SAMPLE_THRESHOLD}" --DP_THRESHOLD "${DP_THRESHOLD}" --ALT_THRESHOLD "${ALT_THRESHOLD}"  --REMOVE_MULTIALLELIC "${REMOVE_MULTIALLELIC}" --PASS "${PASS}" --REMOVE_MITOCHONDRIAL_DNA "${REMOVE_MITOCHONDRIAL_DNA}
python3 "/data/project/Meningioma/script/03.Variant_calling_Annotation/mutect_pipe_hardfilter.py" --INPUT_VCF ${OUTPUT_FMC_PATH} --OUTPUT_VCF ${OUTPUT_FMC_HF_PATH} --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD}  --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA}

#4. Repeat region 지우기
echo "Removing repeat region start: $(date)"
bedtools intersect -header -v -a ${OUTPUT_FMC_HF_PATH} -b ${BLACKLIST} > ${OUTPUT_FMC_HF_RMBLACK_PATH}

#5. STR;STRQ= 있는 line 제거하기 (chr22 제외)
echo "removing STR;STRQ= start: $(date)"
awk '!/STR;STRQ=/ || ($1 == "chr22")' ${OUTPUT_FMC_HF_RMBLACK_PATH} > ${OUTPUT_FMC_HF_RMBLACK_PATH}".temp"
mv ${OUTPUT_FMC_HF_RMBLACK_PATH}".temp" ${OUTPUT_FMC_HF_RMBLACK_PATH}


# bgzip & tabix
echo "bgzip & tabix start: $(date)"
bgzip -c -f ${OUTPUT_FMC_HF_RMBLACK_PATH} > ${OUTPUT_FMC_HF_RMBLACK_PATH}".gz"
tabix -p vcf ${OUTPUT_FMC_HF_RMBLACK_PATH}".gz"
cp ${OUTPUT_FMC_HF_RMBLACK_PATH}".gz.tbi" ${OUTPUT_FMC_HF_RMBLACK_PATH}".tbi"



# VEP
bash "/data/project/Meningioma/script/03.Variant_calling_Annotation/mutect_pair_pipe_20.vep_hg38.sh" \
    --REF ${REF} \
    --INPUT_VCF ${OUTPUT_FMC_HF_RMBLACK_PATH} \
    --OUTPUT_VCF ${OUTPUT_FMC_HF_RMBLACK_PATH%vcf}"vep.vcf"

echo "VEP done: $(date)"





# bgzip -c -f ${OUTPUT_FMC_PATH%".gz"} > ${OUTPUT_FMC_PATH}
# tabix -p vcf ${OUTPUT_FMC_PATH
