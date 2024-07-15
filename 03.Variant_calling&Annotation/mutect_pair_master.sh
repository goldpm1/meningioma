#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

PROJECT_DIR="/data/project/Meningioma"
BAM_DIR=${PROJECT_DIR}"/02.Align"
GVCF_DIR=${PROJECT_DIR}"/05.gvcf"

PON="/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
#REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
REF="/home/goldpm1/reference/genome.fa"
hg="hg38"
gnomad="/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
TMP_PATH=${BAM_DIR}"/temp"
INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi

#for sublog in 01.MTcall 02.FMC_HF_RMBLACK ; do
for sublog in 03.Tumor_pos 04.Other_rescue 05.Bcftools_merge 06.vep; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done
    

sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬


#for idx in ${!sample_name_LIST[@]}; do
#    Sample_ID=${sample_name_LIST[idx]}
for Sample_ID in 190426 220930 221026 221102 221202 230127 230303 230323_2 230405_2 230419 230526 230822 230920; do
    HOLD_J=""
    BCFTOOLS_MERGE_TXT=${BAM_DIR%/*}"/04.mutect/07.2D_merged/01.BCFTOOLS_MERGE_TXT/"${Sample_ID}".txt"
    BCFTOOLS_MERGE_VCF_GZ=${BAM_DIR%/*}"/04.mutect/07.2D_merged/02.BCFTOOLS_MERGE_VCF/"${Sample_ID}".BCFTOOLS_MERGE.vcf.gz"
    BCFTOOLS_MERGE_VCF=${BAM_DIR%/*}"/04.mutect/07.2D_merged/02.BCFTOOLS_MERGE_VCF/"${Sample_ID}".BCFTOOLS_MERGE.vcf"
    BCFTOOLS_MERGE_VEP_VCF_GZ=${BAM_DIR%/*}"/04.mutect/07.2D_merged/03.VEP/"${Sample_ID}".BCFTOOLS_MERGE.vep.vcf.gz"
    BCFTOOLS_MERGE_VEP_VCF=${BAM_DIR%/*}"/04.mutect/07.2D_merged/03.VEP/"${Sample_ID}".BCFTOOLS_MERGE.vep.vcf"
    
    for TISSUE in Tumor Tumor_FFT Tumor_PCT Tumor_PFT Tumor_PP Tumor_PT Dura Ventricle Cortex ; do   #Tumor Dura
        CASE_BAM_PATH=${BAM_DIR}"/"${hg}"/"${TISSUE%_*}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
        CONTROL_BAM_PATH=${BAM_DIR}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
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
            echo $CASE_BAM_PATH

            SAMPLE_THRESHOLD="Dura,Tumor,Cortex,Ventricle"   # "all"
            DP_THRESHOLD=30
            ALT_THRESHOLD=1
            REMOVE_MULTIALLELIC="True"
            PASS="True"
            REMOVE_MITOCHONDRIAL_DNA="True"
            BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
            
            for folder in  ${OUTPUT_VCF_GZ%/*}   ${OUTPUT_FMC_PATH%/*}  ${VEP_FMC_HF_PATH%/*}  ${TUMOR_INTERVAL%/*} ${RESCUE_VCF%/*} ${OTHER_SHARED_VARIANT_VCF%/*} ${TUMOR_UNIQUE_VCF%/*} ${OTHER_UNIQUE_VCF%/*}  ${BCFTOOLS_MERGE_TXT%/*} ${BCFTOOLS_MERGE_VCF%/*} ${BCFTOOLS_MERGE_VEP_VCF_GZ%/*} ${BCFTOOLS_MERGE_VEP_VCF%/*}; do
                if [ ! -d $folder ] ; then
                    mkdir $folder
                fi
            done

            #01. Mutect2 call
            # qsub -pe smp 6 -e $logPath"/01.MTcall" -o $logPath"/01.MTcall" -N 'MT_01.'${Sample_ID}"_"${TISSUE} -hold_jid "doc_"${Sample_ID}"_"${TISSUE} ${CURRENT_PATH}"/mutect_pair_pipe_01.call.sh" \
            # --Sample_ID ${Sample_ID} --CASE_BAM_PATH ${CASE_BAM_PATH} --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} \
            # --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ}  \
            # --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH}

            #02. FMC & HF & RMBLACK & VEP
            # qsub -pe smp 5 -e $logPath"/02.FMC_HF_RMBLACK" -o $logPath"/02.FMC_HF_RMBLACK" -N 'MT_02.'${Sample_ID}"_"${TISSUE} -hold_jid  'MT_01.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/mutect_pair_pipe_02.FMC_HF_RMBLACK.sh" \
            #     --Sample_ID ${Sample_ID} \
            #     --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} \
            #     --OUTPUT_FMC_PATH ${OUTPUT_FMC_PATH} --OUTPUT_FMC_HF_PATH ${OUTPUT_FMC_HF_PATH}  --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
            #     --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH} \
            #     --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} \
            #     --BLACKLIST ${BLACKLIST}

            # 03. Tumor의 경우 Position을 interval로 뽑기  (conda deactivate 필요))   190426은 Tumor_PT를 대상으로 shared mutation 뽑는다 (single cell을 거기서 뽑았기 때문에)
            if [[ "${TISSUE}" == "Tumor" || "${TISSUE}" == "Tumor_PT" ]]; then
                HOLD_J=${HOLD_J}",MT_03."${Sample_ID}"_"${TISSUE}
                qsub -pe smp 1 -e $logPath"/03.Tumor_pos" -o $logPath"/03.Tumor_pos" -N 'MT_03.'${Sample_ID}"_Tumor" -hold_jid  'MT_02.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/mutect_pair_pipe_03.Tumor_position.sh" \
                    --Sample_ID ${Sample_ID} --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
                    --TUMOR_INTERVAL ${TUMOR_INTERVAL} \
                    --BCFTOOLS_MERGE_TXT ${BCFTOOLS_MERGE_TXT}
            elif [[ "${TISSUE}" == "Dura" || "${TISSUE}" == "Ventricle" || "${TISSUE}" == "Cortex"  ]]; then  # 04. (Tumor부터 돌리고) pysam을 기반으로 1개라도 있으면 rescue 해주기
                for folder in  ${RESCUE_VCF%/*} ${OTHER_SHARED_VARIANT_VCF%/*} ${TUMOR_UNIQUE_VCF%/*} ${OTHER_UNIQUE_VCF%/*} ; do
                    if [ -d $folder ] ; then
                        rm -rf $folder
                    fi
                    if [ ! -d $folder ] ; then
                        mkdir $folder
                    fi
                done
                HOLD_J=${HOLD_J}",MT_04."${Sample_ID}"_"${TISSUE}
                if [ "${Sample_ID}" == "190426" ]; then
                    TUMOR_MUTECT2_VCF=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_Tumor_PT.MT2.FMC.HF.RMBLACK.vcf"
                else
                    TUMOR_MUTECT2_VCF=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.vcf"
                fi
                qsub -pe smp 1 -e $logPath"/04.Other_rescue" -o $logPath"/04.Other_rescue" -N 'MT_04.'${Sample_ID}"_"${TISSUE} -hold_jid  'MT_03.'${Sample_ID}"_Tumor"  ${CURRENT_PATH}"/mutect_pair_pipe_04.Other_rescue.sh" \
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
    
    # 05. BCFTOOLS MERGE (같은 날짜끼리 모아주기)
    qsub -pe smp 1 -e $logPath"/05.Bcftools_merge" -o $logPath"/05.Bcftools_merge" -hold_jid ${HOLD_J} -N 'MT_05.'${Sample_ID}  ${CURRENT_PATH}"/mutect_pair_pipe_05.Bcftools_merge.sh" \
        --BCFTOOLS_MERGE_TXT ${BCFTOOLS_MERGE_TXT} \
        --BCFTOOLS_MERGE_VCF_GZ ${BCFTOOLS_MERGE_VCF_GZ} \
        --BCFTOOLS_MERGE_VCF ${BCFTOOLS_MERGE_VCF} \


    # 06. VEP annotation
    qsub -pe smp 6 -e $logPath"/06.vep" -o $logPath"/06.vep" -N "MT_06."${Sample_ID} -hold_jid 'MT_05.'${Sample_ID}  ${CURRENT_PATH}"/mutect_pair_pipe_20.vep.sh" \
        --REF ${REF} \
        --INPUT_VCF ${BCFTOOLS_MERGE_VCF} \
        --OUTPUT_VCF ${BCFTOOLS_MERGE_VEP_VCF}      

done





# if [[ "${Sample_ID}" == "190426"  ]]; then
#     for TISSUE in Tumor_FFT Tumor_PCT Tumor_PFT Tumor_PP Tumor_PT Dura ; do   #Tumor Dura
#         bash  mutect_pipe_23.manualmaf.sh \
#             --INPUT_VCF "/data/project/Meningioma/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf"  \
#             --OUTPUT_MAF "/data/project/Meningioma/04.mutect/08.maf/03.190426/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.maf" \
#             --SELECTED_DB RegBase
#     done
# fi