#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=$(pwd -P)
logPath="${CURRENT_PATH}/log"

PROJECT_DIR="/data/project/Meningioma"
BAM_DIR="${PROJECT_DIR}/02.Align"
GVCF_DIR="${PROJECT_DIR}/05.gvcf"

PON="/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"    
#REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"    # contig: 196
REF="/home/goldpm1/reference/genome.fa"    # contig: 3367
hg="hg38"
gnomad="/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
TMP_PATH="${BAM_DIR}/temp"
INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"

refine_log_dir() {
    local sublog="$1"
    local sublog_path="${logPath}/${sublog}"
    if [ -e "$sublog_path" ] ; then
        rm -rf "$sublog_path"
    fi
    if [ ! -d "$sublog_path" ] ; then
        mkdir -p "$sublog_path"
    fi
}

HOLD_J=','

sample_name_list=$(cat "${CURRENT_PATH%/*}/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # make into array
#for idx in ${!sample_name_LIST[@]}; do
#    Sample_ID=${sample_name_LIST[idx]}


########## For all combinations, run Mutect2 ##############################

# refine_log_dir "02.FMC_HF_RMBLACK"
# refine_log_dir "03.Tumor_pos"

# for Sample_ID in SMS_1st SMS_2nd SMS_3rd SMS_4th SMS_mets; do
#     for Sample_ID in 190426 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT 241127 241127_FT 241127_FTRt 241127_Para \
#         190524 221026 221102 221202 230127 230303_1 230323_2 230323_11 230419 230526 230802_1 230802_2 230802_3 230812 230920 231006 231025 231101 231206 240110 240325 240403_1 240403_2 240403_3 240412_1 240412_2 240412_3 240612 240911 241023 241023_1 \
#         220930_1 220930_2 230405_2 230822 240320 241023_1 \
#         241016 241211 250207 250207_Ant 250425 \
#         241016 241211 250425 250212 250502 250509 \
#         220930_3 220930_4; do

#         for TISSUE in Tumor Dura Dura_Ant DT Falx Falx_1 Falx_2 Ventricle Cortex; do     # Adding underscore for Falx and Dura variants
#             if [[ $TISSUE == *"Falx"* || $TISSUE == *"Falx_1"* || $TISSUE == *"Falx_2"* || $TISSUE == *"Falx_FT"* || $TISSUE == *"Dura"* || $TISSUE == *"DT"* ]]; then
#                 TISSUE_BROAD="Dura"
#             else
#                 TISSUE_BROAD=${TISSUE%%_*}
#             fi

#             CASE_BAM_PATH="${BAM_DIR}/${hg}/${TISSUE_BROAD}/05.Final_bam/${Sample_ID}_${TISSUE}.bam"
#             case "$Sample_ID" in
#                 "190426_FFT"|"190426_PCT"|"190426_PFT"|"190426_PP"|"190426_PT"|"241127_FT"|"241127_FTRt"|"241127_Para"|"241127_1"|"241127_2"|"250207_Ant"|"220930_1"|"220930_2"|"220930_3"|"220930_4"|"230303_1"|"230323_2"|"230323_11"|"230802_1"|"230802_2"|"230802_3"|"240403_1"|"240403_2"|"240403_3"|"240412_1"|"240412_2"|"240412_3"|"241023_1"|"SMS_1st"|"SMS_2nd"|"SMS_3rd"|"SMS_4th"|"SMS_mets")
#                     CONTROL_BAM_PATH="${BAM_DIR}/${hg}/Blood/05.Final_bam/${Sample_ID%%_*}_Blood.bam"
#                     NORMAL_SAMPLE="${Sample_ID%%_*}_Blood"
#                     ;;
#                 "250206_Occ"|"250206_Parie"|"250728_Ant"|"250728_Post")
#                     CONTROL_BAM_PATH="${BAM_DIR}/${hg}/Blood/05.Final_bam/250206_Blood.bam"
#                     NORMAL_SAMPLE="250206_Blood"
#                     INTERVAL="/home/goldpm1/resources/whole.chromosome.GRCh38.bed"
#                     ;;
#                 *)
#                     CONTROL_BAM_PATH="${BAM_DIR}/${hg}/Blood/05.Final_bam/${Sample_ID}_Blood.bam"
#                     NORMAL_SAMPLE="${Sample_ID}_Blood"
#                     ;;
#             esac

#             BCFTOOLS_MERGE_TXT="${BAM_DIR%/*}/04.mutect/07.2D_merged/01.BCFTOOLS_MERGE_TXT/${Sample_ID}.txt"

#             if [ -f "${CASE_BAM_PATH}" ]; then
#                 SAMPLE_THRESHOLD="all"
#                 if [[ "${Sample_ID}" == *"240403"* || "${Sample_ID}" == *"250206_"* || "${Sample_ID}" == *"250728_"* ]]; then
#                     DP_THRESHOLD=14
#                 else
#                     DP_THRESHOLD=30
#                 fi
#                 ALT_THRESHOLD=1
#                 REMOVE_MULTIALLELIC="True"
#                 PASS="True"
#                 REMOVE_MITOCHONDRIAL_DNA="True"
#                 BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
#                 if [[ "${Sample_ID}" == *"220930"* || "${Sample_ID}" == *"221026"* ||  "${Sample_ID}" == *"230127"* || "${Sample_ID}" == *"230323"*  || "${Sample_ID}" == *"230419"* || "${Sample_ID}" == *"231101"* || "${Sample_ID}" == *"240110"* ]]; then
#                     MIN_BQ=20
#                 elif [[ "${Sample_ID}" == *"241016"*  ]]; then
#                     MIN_BQ=30
#                 else
#                     MIN_BQ=25
#                 fi

#                 if [[ "${TISSUE}" == *"Tumor"* ]]; then
#                     if [[ "${Sample_ID}" == *"241016"* || "${Sample_ID}" == *"250212"* || "${Sample_ID}" == *"220930"* || "${Sample_ID}" == *"221026"* || "${Sample_ID}" == *"SMS"* ]]; then
#                         TLOD_THRESHOLD="12.0"
#                     elif [[ "${Sample_ID}" == *"230419"* ]]; then
#                         TLOD_THRESHOLD="10.0"
#                     else
#                         TLOD_THRESHOLD="7.0"
#                     fi
#                 else
#                     TLOD_THRESHOLD="7.0"
#                 fi

#                 OUTPUT_VCF_GZ="${BAM_DIR%/*}/04.mutect/01.raw/${Sample_ID}_${TISSUE}.vcf.gz"
#                 OUTPUT_FMC_PATH="${BAM_DIR%/*}/04.mutect/02.PASS/${Sample_ID}_${TISSUE}.MT2.FMC.vcf"
#                 OUTPUT_FMC_HF_PATH="${BAM_DIR%/*}/04.mutect/02.PASS/${Sample_ID}_${TISSUE}.MT2.FMC.HF.vcf"
#                 OUTPUT_FMC_HF_RMBLACK_PATH="${BAM_DIR%/*}/04.mutect/02.PASS/${Sample_ID}_${TISSUE}.MT2.FMC.HF.RMBLACK.vcf"

#                 # Ensure output directories exist
#                 for folder in \
#                     "${OUTPUT_VCF_GZ%/*}" "${OUTPUT_FMC_PATH%/*}" "${BCFTOOLS_MERGE_TXT%/*}"; do
#                     if [ ! -d "$folder" ]; then
#                         mkdir -p "$folder"
#                     fi
#                 done

#                 # 02. FMC & HF & RMBLACK & VEP
#                 MT02_JOB=$(sbatch \
#                     --job-name="MT_02.${Sample_ID}_${TISSUE}" --output="${logPath}/02.FMC_HF_RMBLACK/${Sample_ID}_${TISSUE}.log" --error="${logPath}/02.FMC_HF_RMBLACK/${Sample_ID}_${TISSUE}.err" \
#                     --cpus-per-task=5 --mem=14GB --time=4:00:00 --partition=cpu --qos nstumor \
#                     --wrap="bash ${CURRENT_PATH}/mutect_pair_pipe_02.FMC_HF_RMBLACK.sh \
#                         --Sample_ID ${Sample_ID} \
#                         --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} \
#                         --OUTPUT_FMC_PATH ${OUTPUT_FMC_PATH} --OUTPUT_FMC_HF_PATH ${OUTPUT_FMC_HF_PATH} --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
#                         --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH} \
#                         --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --TLOD_THRESHOLD ${TLOD_THRESHOLD} --MIN_BQ ${MIN_BQ} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} \
#                         --BLACKLIST ${BLACKLIST}" )
#                 MT02_JOBID=$(echo "${MT02_JOB}" | awk '{print $NF}')

#                 # 03. If Tumor, get position interval
#                 if [[ "${TISSUE}" == "Tumor" ]]; then
#                     TUMOR_INTERVAL="${BAM_DIR%/*}/04.mutect/03.Tumor_interval/${Sample_ID}_Tumor.MT2.FMC.HF.RMBLACK.bed"
#                     MT03_JOB=$(sbatch  --dependency=afterok:${MT02_JOBID} \
#                         --job-name="MT_03.${Sample_ID}_Tumor" \
#                         --output="${logPath}/03.Tumor_pos/${Sample_ID}_Tumor.log" \
#                         --error="${logPath}/03.Tumor_pos/${Sample_ID}_Tumor.err" \
#                         --cpus-per-task=1 --mem=2GB --time=0:10:00 --partition=cpu --qos nstumor \
#                         --wrap="bash ${CURRENT_PATH}/mutect_pair_pipe_03.Tumor_position.sh \
#                             --Sample_ID ${Sample_ID} --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
#                             --TUMOR_INTERVAL ${TUMOR_INTERVAL} \
#                             --BCFTOOLS_MERGE_TXT ${BCFTOOLS_MERGE_TXT}" )
#                     MT03_JOBID=$(echo "${MT03_JOB}" | awk '{print $NF}')
#                 fi

#             fi

#         done
#     done
# done



for Sample_ID in 190426 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT 241127 241127_FT 241127_FTRt 241127_Para \
    190524 221026 221102 221202 230127 230303_1 230323_2 230323_11 230419 230526 230802_1 230802_2 230802_3 230812 230920 231006 231025 231101 231206 240110 240325 240403_1 240403_2 240403_3 240412_1 240412_2 240412_3 240612 240911 241023 241023_1 \
    220930_1 220930_2 230405_2 230822 240320 241023_1 \
    241016 241211 250207 250207_Ant 250425 \
    241016 241211 250425 250212 250502 250509 \
    220930_3 220930_4; do
#for Sample_ID in 230419; do

    BCFTOOLS_MERGE_TXT="${BAM_DIR%/*}/04.mutect/07.2D_merged/01.BCFTOOLS_MERGE_TXT/${Sample_ID}.txt" 
    > "${BCFTOOLS_MERGE_TXT}"  # # "${BCFTOOLS_MERGE_TXT}" 파일의 내용을 비워서(즉, 0바이트로 만들어서) 새로운 내용을 추가할 준비를 한다.
    echo -e "${BAM_DIR%/*}/04.mutect/02.PASS/${Sample_ID}_Tumor.MT2.FMC.HF.RMBLACK.vcf.gz" >> ${BCFTOOLS_MERGE_TXT} # 일단 tumor를 한줄 추가한다


    for TISSUE in Dura Dura_Ant DT Falx Falx_1 Falx_2 Ventricle Cortex; do     # Adding underscore for Falx and Dura variants
        if [[ $TISSUE == *"Falx"* || $TISSUE == *"Falx_1"* || $TISSUE == *"Falx_2"* || $TISSUE == *"Falx_FT"* || $TISSUE == *"Dura"* || $TISSUE == *"DT"* ]]; then
            TISSUE_BROAD="Dura"
        else
            TISSUE_BROAD=${TISSUE%%_*}
        fi


        TUMOR_BAM_PATH="${BAM_DIR}/${hg}/Tumor/05.Final_bam/${Sample_ID}_Tumor.bam"
        RESCUE_BAM_PATH="${BAM_DIR}/${hg}/${TISSUE_BROAD}/05.Final_bam/${Sample_ID}_${TISSUE}.bam"
        
        if [ -f "${RESCUE_BAM_PATH}" ]; then
            case "$Sample_ID" in
                "190426_FFT"|"190426_PCT"|"190426_PFT"|"190426_PP"|"190426_PT"|"241127_FT"|"241127_FTRt"|"241127_Para"|"241127_1"|"241127_2"|"250207_Ant"|"220930_1"|"220930_2"|"220930_3"|"220930_4"|"230303_1"|"230323_2"|"230323_11"|"230802_1"|"230802_2"|"230802_3"|"240403_1"|"240403_2"|"240403_3"|"240412_1"|"240412_2"|"240412_3"|"241023_1"|"SMS_1st"|"SMS_2nd"|"SMS_3rd"|"SMS_4th"|"SMS_mets")
                    BLOOD_BAM_PATH="${BAM_DIR}/${hg}/Blood/05.Final_bam/${Sample_ID%%_*}_Blood.bam"
                    ;;
                "250206_Occ"|"250206_Parie"|"250728_Ant"|"250728_Post")
                    BLOOD_BAM_PATH="${BAM_DIR}/${hg}/Blood/05.Final_bam/250206_Blood.bam"
                    INTERVAL="/home/goldpm1/resources/whole.chromosome.GRCh38.bed"
                    ;;
                *)
                    BLOOD_BAM_PATH="${BAM_DIR}/${hg}/Blood/05.Final_bam/${Sample_ID}_Blood.bam"
                    ;;
            esac

            TUMOR_MUTECT2_VCF="${BAM_DIR%/*}/04.mutect/02.PASS/${Sample_ID}_Tumor.MT2.FMC.HF.RMBLACK.vcf"  # /data/project/Meningioma/04.mutect/02.PASS/190426_FFT_Tumor.MT2.FMC.HF.RMBLACK.vcf
            OUTPUT_FMC_HF_RMBLACK_PATH="${BAM_DIR%/*}/04.mutect/02.PASS/${Sample_ID}_${TISSUE}.MT2.FMC.HF.RMBLACK.vcf"
            OTHER_MUTECT2_VCF="${BAM_DIR%/*}/04.mutect/02.PASS/${Sample_ID}_${TISSUE}.MT2.FMC.HF.RMBLACK.vcf"  # /data/project/Meningioma/04.mutect/02.PASS/190426_Dura.MT2.FMC.HF.RMBLACK.vcf
            RESCUE_VCF="${BAM_DIR%/*}/04.mutect/04.Other_rescue/${Sample_ID}_${TISSUE}.MT2.FMC.HF.RMBLACK.rescue.vcf"  # /data/project/Meningioma/04.mutect/04.Other_rescue/190426_FFT_Dura.MT2.FMC.HF.RMBLACK.rescue.vcf
            TUMOR_SHARED_VARIANT_VCF="${BAM_DIR%/*}/04.mutect/05.Shared_variant/${Sample_ID}_Tumor.MT2.FMC.HF.RMBLACK.shared_variant.vcf"  # /data/project/Meningioma/04.mutect/05.Shared_variant/190426_FFT_Tumor.MT2.FMC.HF.RMBLACK.shared_variant.vcf
            OTHER_SHARED_VARIANT_VCF="${BAM_DIR%/*}/04.mutect/05.Shared_variant/${Sample_ID}_${TISSUE}.MT2.FMC.HF.RMBLACK.shared_variant.vcf"  # /data/project/Meningioma/04.mutect/05.Shared_variant/190426_FFT_Dura.MT2.FMC.HF.RMBLACK.shared_variant.vcf
            TUMOR_UNIQUE_VCF="${BAM_DIR%/*}/04.mutect/06.Unique/${Sample_ID}_${TISSUE}.MT2.FMC.HF.RMBLACK.tumor_unique.vcf"  # /data/project/Meningioma/04.mutect/06.Unique/190426_FFT_Dura.MT2.FMC.HF.RMBLACK.tumor_unique.vcf
            OTHER_UNIQUE_VCF="${BAM_DIR%/*}/04.mutect/06.Unique/${Sample_ID}_${TISSUE}.MT2.FMC.HF.RMBLACK.other_unique.vcf"  # /data/project/Meningioma/04.mutect/06.Unique/190426_FFT_Dura.MT2.FMC.HF.RMBLACK.other_unique.vcf
            # 일단 다 지워주기
            for file in  ${RESCUE_VCF} ${OTHER_SHARED_VARIANT_VCF} ${TUMOR_UNIQUE_VCF} ${OTHER_UNIQUE_VCF} ; do
                find "$(dirname "$file")" -type f -name "$(basename "$file")*" -exec rm -f {} +
            done

            TUMOR_INTERVAL="${BAM_DIR%/*}/04.mutect/03.Tumor_interval/${Sample_ID}_Tumor.MT2.FMC.HF.RMBLACK.bed"
            BCFTOOLS_MERGE_TXT="${BAM_DIR%/*}/04.mutect/07.2D_merged/01.BCFTOOLS_MERGE_TXT/${Sample_ID}.txt"

            THRESHOLD_END=4
            if [ "${Sample_ID}_${TISSUE}" == "250212_Falx" ]; then
                THRESHOLD_BQ=20
            elif [ "${Sample_ID}_${TISSUE}" == "241016_DT" ] || [ "${Sample_ID}_${TISSUE}" == "241211_Dura" ]; then
                THRESHOLD_BQ=11
            else
                THRESHOLD_BQ=30
            fi
            THRESHOLD_MULTIALLELIC=1
            REMOVE_STR=False
            THRESHOLD_TUMOR_CLONAL_VAF=0.0002
            MINIMUM_ALT=1
            REMOVE_CLIP=True


            MT04_JOB=$(sbatch --job-name="MT_04.${Sample_ID}_${TISSUE}" --output="${logPath}/04.Other_rescue/${Sample_ID}_${TISSUE}.log" --error="${logPath}/04.Other_rescue/${Sample_ID}_${TISSUE}.err" \
                --cpus-per-task=6 --mem=12GB --time=4:00:00 --partition=cpu --qos nstumor \
                --wrap="bash ${CURRENT_PATH}/mutect_pair_pipe_04.Other_rescue.sh \
                    --SCRIPT_DIR ${CURRENT_PATH} \
                    --REF ${REF} \
                    --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} --MINIMUM_ALT 1 \
                    --TUMOR_BAM_PATH ${TUMOR_BAM_PATH} \
                    --RESCUE_BAM_PATH ${RESCUE_BAM_PATH} \
                    --BLOOD_BAM_PATH ${BLOOD_BAM_PATH} \
                    --TUMOR_INTERVAL ${TUMOR_INTERVAL} \
                    --TUMOR_MUTECT2_VCF ${TUMOR_MUTECT2_VCF} \
                    --OTHER_MUTECT2_VCF ${OUTPUT_FMC_HF_RMBLACK_PATH} \
                    --RESCUE_VCF ${RESCUE_VCF} \
                    --TUMOR_SHARED_VARIANT_VCF ${TUMOR_SHARED_VARIANT_VCF} \
                    --OTHER_SHARED_VARIANT_VCF ${OTHER_SHARED_VARIANT_VCF} \
                    --TUMOR_UNIQUE_VCF ${TUMOR_UNIQUE_VCF} \
                    --OTHER_UNIQUE_VCF ${OTHER_UNIQUE_VCF} \
                    --TUMOR_INTERVAL ${TUMOR_INTERVAL} \
                    --BCFTOOLS_MERGE_TXT ${BCFTOOLS_MERGE_TXT} \
                    --THRESHOLD_END ${THRESHOLD_END} \
                    --THRESHOLD_BQ ${THRESHOLD_BQ} \
                    --THRESHOLD_MULTIALLELIC ${THRESHOLD_MULTIALLELIC} \
                    --REMOVE_STR ${REMOVE_STR} \
                    --THRESHOLD_TUMOR_CLONAL_VAF ${THRESHOLD_TUMOR_CLONAL_VAF} \
                    --MINIMUM_ALT ${MINIMUM_ALT} \
                    --REMOVE_CLIP ${REMOVE_CLIP}")
            MT04_JOBID=$(echo "${MT04_JOB}" | awk '{print $NF}')
            echo -e "MT_04.${Sample_ID}_${TISSUE} (${MT04_JOBID}) submitted"

        fi
    done

    # 05. BCFTOOLS MERGE (같은 날짜끼리 모아주기)  /data/project/Meningioma/script/03.Variant_calling_Annotation/mutect_pair_pipe_05.Bcftools_merge.sh   #         
    BCFTOOLS_MERGE_VCF="${BAM_DIR%/*}/04.mutect/07.2D_merged/02.BCFTOOLS_MERGE_VCF/${Sample_ID}.BCFTOOLS_MERGE.vcf"
    BCFTOOLS_MERGE_VCF_GZ="${BAM_DIR%/*}/04.mutect/07.2D_merged/02.BCFTOOLS_MERGE_VCF/${Sample_ID}.BCFTOOLS_MERGE.vcf.gz"
    MT05_JOB=$(sbatch --dependency=afterok:${MT04_JOBID} \
        --job-name="MT_05.${Sample_ID}" --output="${logPath}/05.Bcftools_merge/${Sample_ID}.log" --error="${logPath}/05.Bcftools_merge/${Sample_ID}.err" \
        --cpus-per-task=1 --mem=12GB --time=0:10:00 --partition=cpu --qos nstumor \
        --wrap="bash ${CURRENT_PATH}/mutect_pair_pipe_05.Bcftools_merge.sh \
            --BCFTOOLS_MERGE_TXT ${BCFTOOLS_MERGE_TXT} \
            --BCFTOOLS_MERGE_VCF_GZ ${BCFTOOLS_MERGE_VCF_GZ} \
            --BCFTOOLS_MERGE_VCF ${BCFTOOLS_MERGE_VCF}")
    MT05_JOBID=$(echo "${MT05_JOB}" | awk '{print $NF}')
    echo -e "MT_05.${Sample_ID} (${MT05_JOBID}) submitted"


    # 06. VEP annotation
    BCFTOOLS_MERGE_VCF="${BAM_DIR%/*}/04.mutect/07.2D_merged/02.BCFTOOLS_MERGE_VCF/${Sample_ID}.BCFTOOLS_MERGE.vcf"
    BCFTOOLS_MERGE_VEP_VCF="${BAM_DIR%/*}/04.mutect/07.2D_merged/03.VEP/${Sample_ID}.BCFTOOLS_MERGE.vep.vcf"
    MT06_JOB=$(sbatch --dependency=afterok:${MT05_JOBID} \
        --job-name="MT_06.${Sample_ID}" --output="${logPath}/06.vep/${Sample_ID}.log" --error="${logPath}/06.vep/${Sample_ID}.err" \
        --cpus-per-task=6 --mem=12GB --time=4:00:00 --partition=cpu --qos nstumor \
        --wrap="bash ${CURRENT_PATH}/mutect_pair_pipe_20.vep_hg38.sh \
            --REF ${REF} \
            --INPUT_VCF ${BCFTOOLS_MERGE_VCF} \
            --OUTPUT_VCF ${BCFTOOLS_MERGE_VEP_VCF}")
    MT06_JOBID=$(echo "${MT06_JOB}" | awk '{print $NF}')
    echo -e "MT_06.${Sample_ID} (${MT06_JOBID}) submitted"
done