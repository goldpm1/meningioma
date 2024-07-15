#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

PROJECT_DIR="/data/project/Meningioma/61.Lowinput/01.XT_HS"
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

for sublog in 01.MTcall 02.FMC_HF_RMBLACK 03.Tumor_pos 04.Other_rescue 05.Bcftools_merge 06.vep; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


for Sample_ID in 190426; do
    for Clone_No in 1 2; do
        CASE_BAM_PATH=${BAM_DIR}"/"${hg}"/Tumor/06.Final_bam/"${Sample_ID}"_"${Clone_No}".bam"
        CONTROL_BAM_PATH=${BAM_DIR}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
        OUTPUT_VCF_GZ=${BAM_DIR%/*}"/04.mutect/01.raw/"${Sample_ID}"_"${Clone_No}".vcf.gz"
        OUTPUT_FMC_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${Clone_No}".MT2.FMC.vcf"
        OUTPUT_FMC_HF_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${Clone_No}".MT2.FMC.HF.vcf"
        OUTPUT_FMC_HF_RMBLACK_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${Clone_No}".MT2.FMC.HF.RMBLACK.vcf"
        TUMOR_INTERVAL=${BAM_DIR%/*}"/04.mutect/03.Tumor_interval/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.bed"
        # HC_GVCF=${GVCF_DIR}"/02.remove_nonref/"${Sample_ID}"/"${Sample_ID}"_"${Clone_No}".g.vcf"
        # RESCUE_VCF=${BAM_DIR%/*}"/04.mutect/04.Other_rescue/"${Sample_ID}"_"${Clone_No}".MT2.FMC.HF.RMBLACK.rescue.vcf"
        # TUMOR_SHARED_VARIANT_VCF=${BAM_DIR%/*}"/04.mutect/05.Shared_variant/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.shared_variant.vcf"
        # OTHER_SHARED_VARIANT_VCF=${BAM_DIR%/*}"/04.mutect/05.Shared_variant/"${Sample_ID}"_"${Clone_No}".MT2.FMC.HF.RMBLACK.shared_variant.vcf"
        # TUMOR_UNIQUE_VCF=${BAM_DIR%/*}"/04.mutect/06.Unique/"${Sample_ID}"_"${Clone_No}".MT2.FMC.HF.RMBLACK.tumor_unique.vcf"
        # OTHER_UNIQUE_VCF=${BAM_DIR%/*}"/04.mutect/06.Unique/"${Sample_ID}"_"${Clone_No}".MT2.FMC.HF.RMBLACK.other_unique.vcf"

        if [ -f ${CASE_BAM_PATH} ]; then     # File이 있어야만 진행
            echo $CASE_BAM_PATH

            SAMPLE_THRESHOLD="all"
            DP_THRESHOLD=30
            ALT_THRESHOLD=1
            REMOVE_MULTIALLELIC="True"
            PASS="True"
            REMOVE_MITOCHONDRIAL_DNA="True"
            BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
            
            for folder in  ${OUTPUT_VCF_GZ%/*}   ${OUTPUT_FMC_PATH%/*}  ${VEP_FMC_HF_PATH%/*} ${BAM_DIR%/*}"/04.mutect/08.maf/" ${TUMOR_INTERVAL%/*} ${RESCUE_VCF%/*} ${OTHER_SHARED_VARIANT_VCF%/*} ${TUMOR_UNIQUE_VCF%/*} ${OTHER_UNIQUE_VCF%/*}  ${BCFTOOLS_MERGE_TXT%/*} ${BCFTOOLS_MERGE_VCF%/*} ${BCFTOOLS_MERGE_VEP_VCF_GZ%/*} ${BCFTOOLS_MERGE_VEP_VCF%/*}; do
                if [ ! -d $folder ] ; then
                    mkdir -p $folder
                fi
            done

            #01. Mutect2 call
            qsub -pe smp 6 -e $logPath"/01.MTcall" -o $logPath"/01.MTcall" -N 'MT_01.'${Sample_ID}"_"${Clone_No} -hold_jid "doc_"${Sample_ID}"_"${Clone_No} ${CURRENT_PATH}"/03.mutect_pair_pipe_01.call.sh" \
            --Sample_ID ${Sample_ID} --CASE_BAM_PATH ${CASE_BAM_PATH} --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} \
            --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ}  \
            --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH}

            #02. FMC & HF & RMBLACK & VEP
            qsub -pe smp 5 -e $logPath"/02.FMC_HF_RMBLACK" -o $logPath"/02.FMC_HF_RMBLACK" -N 'MT_02.'${Sample_ID}"_"${Clone_No} -hold_jid  'MT_01.'${Sample_ID}"_"${Clone_No}  ${CURRENT_PATH}"/03.mutect_pair_pipe_02.FMC_HF_RMBLACK.sh" \
                --Sample_ID ${Sample_ID} \
                --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} \
                --OUTPUT_FMC_PATH ${OUTPUT_FMC_PATH} --OUTPUT_FMC_HF_PATH ${OUTPUT_FMC_HF_PATH}  --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
                --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH} \
                --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} \
                --BLACKLIST ${BLACKLIST}

        fi
    done
done




    # bash  "mutect_pipe_23.manualmaf.sh" \
    #     --INPUT_VCF ${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${Clone_No}".MT2.FMC.HF.RMBLACK.vep.vcf"  \
    #     --OUTPUT_MAF ${BAM_DIR%/*}"/04.mutect/08.maf/"${Sample_ID}"_"${Clone_No}".MT2.FMC.HF.RMBLACK.vep.maf" \
    #     --SELECTED_DB RegBase

