#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long SCRIPT_DIR:,REF:,Sample_ID:,TISSUE:,TUMOR_BAM_PATH:,RESCUE_BAM_PATH:,BLOOD_BAM_PATH:,TUMOR_INTERVAL:,TUMOR_MUTECT2_VCF:,OTHER_MUTECT2_VCF:,RESCUE_VCF:,TUMOR_SHARED_VARIANT_VCF:,OTHER_SHARED_VARIANT_VCF:,TUMOR_UNIQUE_VCF:,OTHER_UNIQUE_VCF:,BCFTOOLS_MERGE_TXT:,THRESHOLD_END:,THRESHOLD_BQ:,THRESHOLD_MULTIALLELIC:,REMOVE_STR:,THRESHOLD_TUMOR_CLONAL_VAF:,MINIMUM_ALT:,REMOVE_CLIP:, -- "$@")
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
        --REF)
            REF=$2
        shift 2 ;;
        --Sample_ID)
            Sample_ID=$2
        shift 2 ;;
        --TISSUE)
            TISSUE=$2
        shift 2 ;;
        --TUMOR_BAM_PATH)
            TUMOR_BAM_PATH=$2
        shift 2 ;;
        --RESCUE_BAM_PATH)
            RESCUE_BAM_PATH=$2
        shift 2 ;;
        --BLOOD_BAM_PATH)
            BLOOD_BAM_PATH=$2
        shift 2 ;;
        --TUMOR_INTERVAL)
            TUMOR_INTERVAL=$2
        shift 2 ;;
        --TUMOR_MUTECT2_VCF)
            TUMOR_MUTECT2_VCF=$2
        shift 2 ;;
        --OTHER_MUTECT2_VCF)
            OTHER_MUTECT2_VCF=$2
        shift 2 ;;
        --RESCUE_VCF)
            RESCUE_VCF=$2
        shift 2 ;;
        --TUMOR_SHARED_VARIANT_VCF)
            TUMOR_SHARED_VARIANT_VCF=$2
        shift 2 ;;
        --OTHER_SHARED_VARIANT_VCF)
            OTHER_SHARED_VARIANT_VCF=$2
        shift 2 ;;
        --TUMOR_UNIQUE_VCF)
            TUMOR_UNIQUE_VCF=$2
        shift 2 ;;
        --OTHER_UNIQUE_VCF)
            OTHER_UNIQUE_VCF=$2
        shift 2 ;;
        --BCFTOOLS_MERGE_TXT)
            BCFTOOLS_MERGE_TXT=$2
        shift 2 ;;
        --THRESHOLD_END)
            THRESHOLD_END=$2
        shift 2 ;;
        --THRESHOLD_BQ)
            THRESHOLD_BQ=$2
        shift 2 ;;
        --THRESHOLD_MULTIALLELIC)
            THRESHOLD_MULTIALLELIC=$2
        shift 2 ;;
        --REMOVE_STR)
            REMOVE_STR=$2
        shift 2 ;;
        --THRESHOLD_TUMOR_CLONAL_VAF)
            THRESHOLD_TUMOR_CLONAL_VAF=$2
        shift 2 ;;
        --MINIMUM_ALT)
            MINIMUM_ALT=$2
        shift 2 ;;
        --REMOVE_CLIP)
            REMOVE_CLIP=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo -e "Start time : "$(date) "\n\n"

echo -e "Sample_ID : "${Sample_ID}
echo -e "TISSUE : "${TISSUE}
echo -e "TUMOR_BAM_PATH : "${TUMOR_BAM_PATH}
echo -e "RESCUE_BAM_PATH : "${RESCUE_BAM_PATH}
echo -e "BLOOD_BAM_PATH : "${BLOOD_BAM_PATH}
echo -e "TUMOR_INTERVAL : "${TUMOR_INTERVAL}
echo -e "TUMOR_MUTECT2_VCF : "${TUMOR_MUTECT2_VCF}
echo -e "OTHER_MUTECT2_VCF : "${OTHER_MUTECT2_VCF}
echo -e "RESCUE_VCF : "${RESCUE_VCF}
echo -e "TUMOR_SHARED_VARIANT_VCF : "${TUMOR_SHARED_VARIANT_VCF}
echo -e "OTHER_SHARED_VARIANT_VCF : "${OTHER_SHARED_VARIANT_VCF}
echo -e "TUMOR_UNIQUE_VCF : "${TUMOR_UNIQUE_VCF}
echo -e "OTHER_UNIQUE_VCF : "${OTHER_UNIQUE_VCF}
echo -e "BCFTOOLS_MERGE_TXT : "${BCFTOOLS_MERGE_TXT}
echo -e "THRESHOLD_END : "${THRESHOLD_END}
echo -e "THRESHOLD_BQ : "${THRESHOLD_BQ}
echo -e "THRESHOLD_MULTIALLELIC : "${THRESHOLD_MULTIALLELIC}
echo -e "REMOVE_STR : "${REMOVE_STR}
echo -e "THRESHOLD_TUMOR_CLONAL_VAF : "${THRESHOLD_TUMOR_CLONAL_VAF}
echo -e "MINIMUM_ALT : "${MINIMUM_ALT}
echo -e "REMOVE_CLIP : "${REMOVE_CLIP}


# 일단 복사해놓고 그 뒤에 rescue 를 추가한다
cp -rf ${OTHER_MUTECT2_VCF}  ${RESCUE_VCF}


# header 를 복사해놓기
#echo -e "\n\ngrep '#' ${TUMOR_MUTECT2_VCF} > ${TUMOR_SHARED_VARIANT_VCF}"
grep '#' ${TUMOR_MUTECT2_VCF} > ${TUMOR_SHARED_VARIANT_VCF}
grep '#' ${OTHER_MUTECT2_VCF} > ${OTHER_SHARED_VARIANT_VCF}
grep '#' ${TUMOR_MUTECT2_VCF} > ${TUMOR_UNIQUE_VCF}
grep '#' ${OTHER_MUTECT2_VCF} > ${OTHER_UNIQUE_VCF}

echo -e "python3 ${SCRIPT_DIR}"/mutect_pair_pipe_04.Other_rescue.py" \
    --Sample_ID ${Sample_ID} \
    --TISSUE ${TISSUE} \
    --TUMOR_INTERVAL ${TUMOR_INTERVAL} \
    --TUMOR_BAM_PATH ${TUMOR_BAM_PATH} \
    --RESCUE_BAM_PATH ${RESCUE_BAM_PATH} \
    --BLOOD_BAM_PATH ${BLOOD_BAM_PATH} \
    --TUMOR_MUTECT2_VCF ${TUMOR_MUTECT2_VCF} \
    --OTHER_MUTECT2_VCF ${OTHER_MUTECT2_VCF} \
    --RESCUE_VCF ${RESCUE_VCF} \
    --TUMOR_SHARED_VARIANT_VCF ${TUMOR_SHARED_VARIANT_VCF} \
    --OTHER_SHARED_VARIANT_VCF ${OTHER_SHARED_VARIANT_VCF} \
    --TUMOR_UNIQUE_VCF ${TUMOR_UNIQUE_VCF} \
    --OTHER_UNIQUE_VCF ${OTHER_UNIQUE_VCF} \
    --THRESHOLD_END ${THRESHOLD_END} \
    --THRESHOLD_BQ ${THRESHOLD_BQ} \
    --THRESHOLD_MULTIALLELIC ${THRESHOLD_MULTIALLELIC} \
    --REMOVE_STR ${REMOVE_STR} \
    --THRESHOLD_TUMOR_CLONAL_VAF ${THRESHOLD_TUMOR_CLONAL_VAF} \
    --MINIMUM_ALT ${MINIMUM_ALT} \
    --REMOVE_CLIP ${REMOVE_CLIP}
"


python3 ${SCRIPT_DIR}"/mutect_pair_pipe_04.Other_rescue.py" \
    --Sample_ID ${Sample_ID} \
    --TISSUE ${TISSUE} \
    --TUMOR_INTERVAL ${TUMOR_INTERVAL} \
    --TUMOR_BAM_PATH ${TUMOR_BAM_PATH} \
    --RESCUE_BAM_PATH ${RESCUE_BAM_PATH} \
    --BLOOD_BAM_PATH ${BLOOD_BAM_PATH} \
    --TUMOR_MUTECT2_VCF ${TUMOR_MUTECT2_VCF} \
    --OTHER_MUTECT2_VCF ${OTHER_MUTECT2_VCF} \
    --RESCUE_VCF ${RESCUE_VCF} \
    --TUMOR_SHARED_VARIANT_VCF ${TUMOR_SHARED_VARIANT_VCF} \
    --OTHER_SHARED_VARIANT_VCF ${OTHER_SHARED_VARIANT_VCF} \
    --TUMOR_UNIQUE_VCF ${TUMOR_UNIQUE_VCF} \
    --OTHER_UNIQUE_VCF ${OTHER_UNIQUE_VCF} \
    --THRESHOLD_END ${THRESHOLD_END} \
    --THRESHOLD_BQ ${THRESHOLD_BQ} \
    --THRESHOLD_MULTIALLELIC ${THRESHOLD_MULTIALLELIC} \
    --REMOVE_STR ${REMOVE_STR} \
    --THRESHOLD_TUMOR_CLONAL_VAF ${THRESHOLD_TUMOR_CLONAL_VAF} \
    --MINIMUM_ALT ${MINIMUM_ALT} \
    --REMOVE_CLIP ${REMOVE_CLIP}


# sort lexicographically
grep "^#" ${RESCUE_VCF} > ${RESCUE_VCF}".sorted" && grep -v "^#" ${RESCUE_VCF} | sort -k1,1V -k2n >> ${RESCUE_VCF}".sorted"
mv ${RESCUE_VCF}".sorted" ${RESCUE_VCF}


# Subtract variants in file2.vcf from file1.vcf
vcftools --vcf ${TUMOR_MUTECT2_VCF} --exclude-positions ${OTHER_SHARED_VARIANT_VCF}  --recode --out ${TUMOR_UNIQUE_VCF}
vcftools --vcf ${RESCUE_VCF} --exclude-positions ${OTHER_SHARED_VARIANT_VCF}  --recode --out ${OTHER_UNIQUE_VCF}

mv ${TUMOR_UNIQUE_VCF}".recode.vcf" ${TUMOR_UNIQUE_VCF}
rm -rf ${TUMOR_UNIQUE_VCF}".recode.vcf"
mv ${OTHER_UNIQUE_VCF}".recode.vcf" ${OTHER_UNIQUE_VCF}
rm -rf ${OTHER_UNIQUE_VCF}".recode.vcf"



# merge shared variant
# bgzip & tabix
bgzip -c -f ${RESCUE_VCF} > ${RESCUE_VCF}".gz"
tabix -p vcf ${RESCUE_VCF}".gz"
bgzip -c -f ${TUMOR_SHARED_VARIANT_VCF} > ${TUMOR_SHARED_VARIANT_VCF}".gz"
tabix -p vcf ${TUMOR_SHARED_VARIANT_VCF}".gz"
bgzip -c -f ${OTHER_SHARED_VARIANT_VCF} > ${OTHER_SHARED_VARIANT_VCF}".gz"
tabix -p vcf ${OTHER_SHARED_VARIANT_VCF}".gz"
bgzip -c -f ${TUMOR_UNIQUE_VCF} > ${TUMOR_UNIQUE_VCF}".gz"
tabix -p vcf ${TUMOR_UNIQUE_VCF}".gz"
bgzip -c -f ${OTHER_UNIQUE_VCF} > ${OTHER_UNIQUE_VCF}".gz"
tabix -p vcf ${OTHER_UNIQUE_VCF}".gz"


# VEP
#for VCF in ${RESCUE_VCF} ${TUMOR_SHARED_VARIANT_VCF} ${OTHER_SHARED_VARIANT_VCF} ${TUMOR_UNIQUE_VCF} ${OTHER_UNIQUE_VCF}; do
for VCF in ${RESCUE_VCF} ; do
    bash "/data/project/Meningioma/script/03.Variant_calling_Annotation/mutect_pair_pipe_20.vep_hg38.sh" \
        --REF ${REF} \
        --INPUT_VCF ${VCF} \
        --OUTPUT_VCF ${VCF%vcf}"vep.vcf"
done


# bcftools merge txt 파일에 덧붙이기
echo -e "echo ${RESCUE_VCF}".gz" >> ${BCFTOOLS_MERGE_TXT}"
echo ${RESCUE_VCF}".gz" >> ${BCFTOOLS_MERGE_TXT}



echo -e "\n\nEnd time : "$(date)