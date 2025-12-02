#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

PROJECT_DIR="/data/project/Meningioma"
BAM_DIR=${PROJECT_DIR}"/02.Align"
GVCF_DIR=${PROJECT_DIR}"/05.gvcf"

PON="/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"    
#REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"    # contig: 196
REF="/home/goldpm1/reference/genome.fa"    # contig: 3367
hg="hg38"
gnomad="/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
TMP_PATH=${BAM_DIR}"/temp"
INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi

for sublog in 01.MTcall 02.FMC_HF_RMBLACK 03.Tumor_pos 04.Other_rescue 05.Bcftools_merge 06.vep; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done
    
HOLD_J=','
sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬
#for idx in ${!sample_name_LIST[@]}; do
#    Sample_ID=${sample_name_LIST[idx]}


################### 읻단 모든 조합에 대해 Mutect2 돌리기 ##############################33

for Sample_ID in 190426 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT 241127 241127_FT 241127_FTRt 241127_Para \
    190524 221026 221102 221202 230127 230303_1 230323_2 230323_11 230419 230526 230802_1 230802_2 230802_3 230812 230920 231006 231025 231101 231206  240110 240325 240403_1 240403_2 240403_3 240412_1 240412_2 240412_3 240612 240911 241023 241023_1 ; do     
      
    for TISSUE in Tumor Dura Falx_1 Falx_2 Ventricle Cortex ; do            # Dura, Falx는 여기에다가 _를 넣음 (Falx_1, Falx_2)
        if [[ $TISSUE == *"Falx_1"* || $TISSUE == *"Falx_2"* ]]; then
            TISSUE_BROAD="Dura"
        else
            TISSUE_BROAD=${TISSUE%%_*}
        fi

        CASE_BAM_PATH=${BAM_DIR}"/"${hg}"/"${TISSUE_BROAD}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"       # 05.Final_bam  vs  #01.Pre_bam/ *.sorted.bam
        
        case "$Sample_ID" in
            "190426_FFT"|"190426_PCT"|"190426_PFT"|"190426_PP"|"190426_PT"|"241023_1"|"241127_FT"|"241127_FTRt"|"241127_Para"|"241127_1"|"241127_2")
                CONTROL_BAM_PATH=${BAM_DIR}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID%%_*}"_Blood.bam"
                NORMAL_SAMPLE=${Sample_ID%%_*}"_Blood"
                ;;
            *)
                
                CONTROL_BAM_PATH=${BAM_DIR}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
                NORMAL_SAMPLE=${Sample_ID}"_Blood"
                ;;
            esac
        OUTPUT_VCF_GZ=${BAM_DIR%/*}"/04.mutect/01.raw/"${Sample_ID}"_"${TISSUE}".vcf.gz"
        OUTPUT_FMC_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.vcf"
        OUTPUT_FMC_HF_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.vcf"
        OUTPUT_FMC_HF_RMBLACK_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vcf"
        TUMOR_INTERVAL=${BAM_DIR%/*}"/04.mutect/03.Tumor_interval/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.bed"
        HC_GVCF=${GVCF_DIR}"/02.remove_nonref/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".g.vcf"
        RESCUE_VCF=${BAM_DIR%/*}"/04.mutect/04.Other_rescue/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vcf"
        TUMOR_SHARED_VARIANT_VCF=${BAM_DIR%/*}"/04.mutect/05.Shared_variant/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.shared_variant.vcf"
        OTHER_SHARED_VARIANT_VCF=${BAM_DIR%/*}"/04.mutect/05.Shared_variant/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.shared_variant.vcf"
        TUMOR_UNIQUE_VCF=${BAM_DIR%/*}"/04.mutect/06.Unique/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.tumor_unique.vcf"
        OTHER_UNIQUE_VCF=${BAM_DIR%/*}"/04.mutect/06.Unique/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.other_unique.vcf"


        if [ -f ${CASE_BAM_PATH} ]; then     # File이 있어야만 진행
            #echo ${CASE_BAM_PATH}

            SAMPLE_THRESHOLD="all"
            if [[ "${Sample_ID}" == *"240403"* ]]; then
                DP_THRESHOLD=14
            else
                DP_THRESHOLD=30
            fi
            ALT_THRESHOLD=1
            REMOVE_MULTIALLELIC="True"
            PASS="True"
            REMOVE_MITOCHONDRIAL_DNA="True"
            BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
            if [[ "${Sample_ID}" == *"220930"* ||  "${Sample_ID}" == *"221026"* ||  "${Sample_ID}" == *"230127"* || "${Sample_ID}" == *"230323"*  || "${Sample_ID}" == *"230419"* || "${Sample_ID}" == *"231101"* || "${Sample_ID}" == *"240110"*  ]]; then
                MIN_BQ=20
            else
                MIN_BQ=25
            fi

            # 디렉토리 없으면 만들어주기            
            for folder in  ${OUTPUT_VCF_GZ%/*}   ${OUTPUT_FMC_PATH%/*}  ${VEP_FMC_HF_PATH%/*}  ${TUMOR_INTERVAL%/*} ${RESCUE_VCF%/*} ${OTHER_SHARED_VARIANT_VCF%/*} ${TUMOR_UNIQUE_VCF%/*} ${OTHER_UNIQUE_VCF%/*}  ${BCFTOOLS_MERGE_TXT%/*} ${BCFTOOLS_MERGE_VCF%/*} ${BCFTOOLS_MERGE_VEP_VCF_GZ%/*} ${BCFTOOLS_MERGE_VEP_VCF%/*}; do
                if [ ! -d $folder ] ; then
                    mkdir $folder
                fi
            done

            # #01. Mutect2 call
            # qsub -pe smp 6 -e $logPath"/01.MTcall" -o $logPath"/01.MTcall" -N 'MT_01.'${Sample_ID}"_"${TISSUE} -hold_jid "postbwa_"${Sample_ID}"_"${TISSUE} ${CURRENT_PATH}"/mutect_pair_pipe_01.call.sh" \
            # --Sample_ID ${Sample_ID} --CASE_BAM_PATH ${CASE_BAM_PATH} --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} \
            # --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ}  \
            # --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --NORMAL_SAMPLE ${NORMAL_SAMPLE} --TMP_PATH ${TMP_PATH}

            #02. FMC & HF & RMBLACK & VEP
            qsub -pe smp 5 -e $logPath"/02.FMC_HF_RMBLACK" -o $logPath"/02.FMC_HF_RMBLACK" -N 'MT_02.'${Sample_ID}"_"${TISSUE} -hold_jid  'MT_01.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/mutect_pair_pipe_02.FMC_HF_RMBLACK.sh" \
                --Sample_ID ${Sample_ID} --TISSUE_BROAD ${TISSUE_BROAD} \
                --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} \
                --OUTPUT_FMC_PATH ${OUTPUT_FMC_PATH} --OUTPUT_FMC_HF_PATH ${OUTPUT_FMC_HF_PATH}  --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
                --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH} \
                --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --MIN_BQ ${MIN_BQ} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} \
                --BLACKLIST ${BLACKLIST}


            # 03. Tumor의 경우 Position을 interval로 뽑기  (conda deactivate 필요)
            if [[ "${TISSUE}" == "Tumor"  ]]; then 
                TUMOR_INTERVAL=${BAM_DIR%/*}"/04.mutect/03.Tumor_interval/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.bed"
                
                qsub -pe smp 1 -e $logPath"/03.Tumor_pos" -o $logPath"/03.Tumor_pos" -N 'MT_03.'${Sample_ID}"_Tumor" -hold_jid  'MT_02.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/mutect_pair_pipe_03.Tumor_position.sh" \
                    --Sample_ID ${Sample_ID} --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
                    --TUMOR_INTERVAL ${TUMOR_INTERVAL} \
                    --BCFTOOLS_MERGE_TXT ${BCFTOOLS_MERGE_TXT}

            elif [[ "${TISSUE_BROAD}" == "Dura" || "${TISSUE}" == "Falx_1" || "${TISSUE}" == "Falx_2" || "${TISSUE}" == "Ventricle" || "${TISSUE}" == "Cortex"   ]]; then  # 04. (Tumor부터 돌리고) pysam을 기반으로 1개라도 있으면 rescue 해주기
                HOLD_J=${HOLD_J}",MT_04."${Sample_ID}"_"${TISSUE}

                # 일단 다 지워주기 (이 명령어는 $file 변수에 저장된 파일의 디렉토리에서, $file 이름으로 시작하는 모든 파일(예: 확장자가 다르거나 뒤에 추가 문자열이 붙은 파일 등)을 찾아서 삭제(rm -f)하는 명령어입니다.)
                for file in  ${RESCUE_VCF} ${OTHER_SHARED_VARIANT_VCF} ${TUMOR_UNIQUE_VCF} ${OTHER_UNIQUE_VCF} ; do
                    find "$(dirname "$file")" -type f -name "$(basename "$file")*" -exec rm -f {} +
                done

                # Matched tumor 지정해주기
                if [ "${Sample_ID}" == "190426" ]; then                # 190426_Dura
                    TUMOR_MUTECT2_VCF=${BAM_DIR%/*}"/04.mutect/02.PASS/190426_Tumor_PT.MT2.FMC.HF.RMBLACK.vcf"
                    HOLD_JJ='MT_02.'${Sample_ID}"_"${TISSUE}",MT_03.190426_PT_Tumor"
                elif [ "${Sample_ID}" == "241127" ]; then    # 241127_Falx_1, 241127_Falx_2              # "${Sample_ID%%_*}" == "241127"
                    TUMOR_MUTECT2_VCF=${BAM_DIR%/*}"/04.mutect/02.PASS/241127_Tumor.MT2.FMC.HF.RMBLACK.vcf"
                    HOLD_JJ='MT_02.'${Sample_ID}"_"${TISSUE}",MT_03.241127_Tumor"
                elif [ -n "${Sample_ID}" ]; then
                    TUMOR_MUTECT2_VCF=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.vcf"
                    HOLD_JJ='MT_02.'${Sample_ID}"_"${TISSUE}',MT_03.'${Sample_ID}"_Tumor"
                else
                    echo "Error: Sample_ID is not set or invalid."
                    exit 1
                fi

                qsub -pe smp 6 -e $logPath"/04.Other_rescue" -o $logPath"/04.Other_rescue" -N 'MT_04.'${Sample_ID}"_"${TISSUE} -hold_jid  ${HOLD_JJ} ${CURRENT_PATH}"/mutect_pair_pipe_04.Other_rescue.sh" \
                    --SCRIPT_DIR ${CURRENT_PATH} \
                    --REF ${REF} \
                    --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} --MINIMUM_ALT 1 \
                    --TUMOR_INTERVAL ${TUMOR_INTERVAL} \
                    --CASE_BAM_PATH ${CASE_BAM_PATH} \
                    --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} \
                    --TUMOR_MUTECT2_VCF ${TUMOR_MUTECT2_VCF} \
                    --OTHER_MUTECT2_VCF ${OUTPUT_FMC_HF_RMBLACK_PATH} \
                    --HC_GVCF ${HC_GVCF} \
                    --RESCUE_VCF ${RESCUE_VCF} \
                    --TUMOR_SHARED_VARIANT_VCF ${TUMOR_SHARED_VARIANT_VCF} \
                    --OTHER_SHARED_VARIANT_VCF ${OTHER_SHARED_VARIANT_VCF} \
                    --TUMOR_UNIQUE_VCF ${TUMOR_UNIQUE_VCF} \
                    --OTHER_UNIQUE_VCF ${OTHER_UNIQUE_VCF} \
                    --BCFTOOLS_MERGE_TXT ${BCFTOOLS_MERGE_TXT}
            fi        
        fi
    done

    HOLD_J="${HOLD_J:1}"  # 맨 앞 ,를 빼줌

    # # 05. BCFTOOLS MERGE (같은 날짜끼리 모아주기)  /data/project/Meningioma/script/03.Variant_calling_Annotation/mutect_pair_pipe_05.Bcftools_merge.sh   # 
    # qsub -pe smp 1 -e $logPath"/05.Bcftools_merge" -o $logPath"/05.Bcftools_merge"  -hold_jid ${HOLD_J}   -N 'MT_05.'${Sample_ID}  ${CURRENT_PATH}"/mutect_pair_pipe_05.Bcftools_merge.sh" \
    #     --BCFTOOLS_MERGE_TXT ${BCFTOOLS_MERGE_TXT} \
    #     --BCFTOOLS_MERGE_VCF_GZ ${BCFTOOLS_MERGE_VCF_GZ} \
    #     --BCFTOOLS_MERGE_VCF ${BCFTOOLS_MERGE_VCF} \


    # # 06. VEP annotation
    # qsub -pe smp 6 -e $logPath"/06.vep" -o $logPath"/06.vep" -N "MT_06."${Sample_ID} -hold_jid 'MT_05.'${Sample_ID}  ${CURRENT_PATH}"/mutect_pair_pipe_20.vep.sh" \
    #     --REF ${REF} \
    #     --INPUT_VCF ${BCFTOOLS_MERGE_VCF} \
    #     --OUTPUT_VCF ${BCFTOOLS_MERGE_VEP_VCF}      

done







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


# ######## 241127 #############
# for Sample_ID in 241127 241127_FT 241127_FTRt 241127_Para; do     
#     for TISSUE in Tumor Falx_1 Falx_2 ; do
#         if [ -f "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf" ]; then                 # 해당 파일이 존재하면
#             bash  mutect_pipe_23.manualmaf.sh \
#                 --INPUT_VCF "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf"  \
#                 --OUTPUT_MAF "/data/project/Meningioma/04.mutect/08.maf/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.maf" \
#                 --SELECTED_DB RegBase
#         fi
#     done
# done
