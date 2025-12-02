#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

PROJECT_DIR="/data/project/Meningioma"
BAM_DIR=${PROJECT_DIR}"/02.Align"
MUTECT_DIR=${PROJECT_DIR}"/04.mutect/04.Other_rescue"




# ####################################### Make MAF for meningiomatosis (여기서는 RESCUE 이전의 VCF를 가지고 수행) ##############################

# for Sample_ID in 190426 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT; do     
#     for TISSUE in Dura Tumor  ; do   #Tumor Dura
#         if [ -f "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf" ]; then                 # 해당 파일이 존재하면
#             bash  mutect_pipe_23.manualmaf.sh \
#                 --INPUT_VCF "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf"  \
#                 --OUTPUT_MAF "/data/project/Meningioma/04.mutect/08.maf/190426/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.maf" \
#                 --SELECTED_DB RegBase
#             fi
#     done
# fi

# ######## 230323 #############
# for Sample_ID in 230323_2 230323_11; do
#     for TISSUE in Tumor Dura ; do
#         bash  mutect_pipe_23.manualmaf.sh \
#             --INPUT_VCF "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf"  \
#             --OUTPUT_MAF "/data/project/Meningioma/04.mutect/08.maf/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.maf" \
#             --SELECTED_DB RegBase
#     done
# done


######## SMS #############
Sample_ID_Representative="SMS"
for Sample_ID in SMS_1st SMS_2nd SMS_3rd SMS_4th SMS_mets; do
    for TISSUE in Tumor ; do
        bash  mutect_pipe_23.manualmaf.sh \
            --INPUT_VCF "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf"  \
            --OUTPUT_MAF "/data/project/Meningioma/04.mutect/08.maf/"${Sample_ID_Representative}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.maf" \
            --SELECTED_DB "RegBase"
        echo -e ${Sample_ID}"_"${TISSUE}" done : "${date}
    done
done




######## 241127 #############
# for Sample_ID in 241127 241127_FT 241127_FTRt 241127_Para; do     
#     for TISSUE in Tumor Falx_1 Falx_2 ; do
#         if [ -f "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf" ]; then                 # 해당 파일이 존재하면
#             echo -e "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf"
#             bash  mutect_pipe_23.manualmaf.sh \
#                 --INPUT_VCF "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf"  \
#                 --OUTPUT_MAF "/data/project/Meningioma/04.mutect/08.maf/"${Sample_ID%%_*}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.maf" \
#                 --SELECTED_DB RegBase
#         fi
#     done
# done

# for Sample_ID in 241127 241127_FT 241127_FTRt 241127_Para; do     
#     for TISSUE in Tumor Falx_1 Falx_2 ; do
#         if [ -f "/data/project/Meningioma/04.mutect/04.Other_rescue/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.vcf" ]; then                 # 해당 파일이 존재하면
#             echo -e "/data/project/Meningioma/04.mutect/04.Other_rescue/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.vcf"
#             bash  mutect_pipe_23.manualmaf.sh \
#                 --INPUT_VCF "/data/project/Meningioma/04.mutect/04.Other_rescue/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.vcf"  \
#                 --OUTPUT_MAF "/data/project/Meningioma/04.mutect/08.maf/"${Sample_ID%%_*}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.maf" \
#                 --SELECTED_DB RegBase
#         fi
#     done
# done


############################### 4개의 rescue.vep.VCF 파일을 하나로 합치고, 중복(동일 chromosome, position)은 1개만 남기고, 정렬까지 수행하는 코드 #######################

# # 파일 경로 패턴 지정
# Sample_ID="241127"
# TISSUE="Falx_1"
# PATTERN="${Sample_ID}*${TISSUE}*.rescue.vep.vcf"
# VCF_LIST=($(ls ${MUTECT_DIR}/${PATTERN} 2>/dev/null))
# #echo "VCF_LIST: ${VCF_LIST[@]}"
# OUTPUT_MERGED_VCF="${MUTECT_DIR}/${Sample_ID}_${TISSUE}.rescue.merged.sorted.vep.vcf"


# # header 추출 (맨 앞 파일에서 header만 추출)
# grep "^#" "${VCF_LIST[0]}" > "${OUTPUT_MERGED_VCF}.tmp"
# # body(variant) 부분 합치기
# for vcf in "${VCF_LIST[@]}"; do
#     grep -v "^#" "$vcf"
# done | awk '!seen[$1,$2]++' >> "${OUTPUT_MERGED_VCF}.tmp"  # awk '!seen[$1,$2]++'는 chromosome($1)과 position($2)이 같은 경우 중복을 제거하고 한 번만 남깁니다.

# # chromosome, position 기준으로 정렬 (header 제외하고 정렬)
# grep "^#" "${OUTPUT_MERGED_VCF}.tmp" > "${OUTPUT_MERGED_VCF}"
# grep -v "^#" "${OUTPUT_MERGED_VCF}.tmp" | sort -k1,1V -k2,2n >> "${OUTPUT_MERGED_VCF}"
# rm -f "${OUTPUT_MERGED_VCF}.tmp"

# # Make MAF for Falx
# echo -e "/data/project/Meningioma/04.mutect/04.Other_rescue/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.vcf"
# bash  mutect_pipe_23.manualmaf.sh \
#     --INPUT_VCF ${OUTPUT_MERGED_VCF}  \
#     --OUTPUT_MAF "/data/project/Meningioma/04.mutect/08.maf/"${Sample_ID%%_*}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.maf" \
#     --SELECTED_DB RegBase

# # Make MAF for each tumour
# for Sample_ID in 241127 241127_FT 241127_FTRt 241127_Para; do     
#     TISSUE="Tumor"
#     if [ -f "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf" ]; then                 # 해당 파일이 존재하면
#         echo -e "/data/project/Meningioma/04.mutect/04.Other_rescue/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.vcf"
#         bash  mutect_pipe_23.manualmaf.sh \
#             --INPUT_VCF "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf"  \
#             --OUTPUT_MAF "/data/project/Meningioma/04.mutect/08.maf/"${Sample_ID%%_*}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.maf" \
#             --SELECTED_DB RegBase
#     fi
# done