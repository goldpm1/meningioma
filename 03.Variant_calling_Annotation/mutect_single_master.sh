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
INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"      # SureSelect V5
INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_V8_hg38/S33613367_Covered.bed"      # SureSelect V8


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
    

#for Sample_ID in 250310 250319 250407 250408 250428 250509 250513 250515 250520 250526_Left 250526_Right 250527 250530  250602 250605 250609  ; do 
for Sample_ID in 220930_3 221021 231011 231208 231220 240112 250227 250617 250627 250704 250716 250723_Ant 250723_Post 250725; do 
    HOLD_J=""
    
    for TISSUE in Tumor ; do   #Tumor Dura
        CASE_BAM_PATH=${BAM_DIR}"/"${hg}"/"${TISSUE%_*}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"       # 05.Final_bam  vs  #01.Pre_bam/ *.sorted.bam
        CONTROL_BAM_PATH=${BAM_DIR}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
        OUTPUT_VCF_GZ=${BAM_DIR%/*}"/04.mutect/01.raw/"${Sample_ID}"_"${TISSUE}".single.vcf.gz"
        OUTPUT_FMC_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".single.MT2.FMC.vcf"
        OUTPUT_FMC_HF_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".single.MT2.FMC.HF.vcf"
        OUTPUT_FMC_HF_RMBLACK_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".single.MT2.FMC.HF.RMBLACK.vcf"


        if [ -f ${CASE_BAM_PATH} ]; then     # File이 있어야만 진행

            SAMPLE_THRESHOLD="Dura,Tumor,Cortex,Ventricle"   # "all"
            DP_THRESHOLD=30
            ALT_THRESHOLD=1
            MIN_BQ=20
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
            qsub -pe smp 6 -e $logPath"/01.MTcall" -o $logPath"/01.MTcall" -N 'MT_01.'${Sample_ID}"_"${TISSUE} -hold_jid "doc_"${Sample_ID}"_"${TISSUE} ${CURRENT_PATH}"/mutect_single_pipe_01.call.sh" \
            --Sample_ID ${Sample_ID} --CASE_BAM_PATH ${CASE_BAM_PATH}  \
            --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ}  \
            --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH}

            #02. FMC & HF & RMBLACK & VEP
            qsub -pe smp 5 -e $logPath"/02.FMC_HF_RMBLACK" -o $logPath"/02.FMC_HF_RMBLACK" -N 'MT_02.'${Sample_ID}"_"${TISSUE} -hold_jid  'MT_01.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/mutect_pair_pipe_02.FMC_HF_RMBLACK.sh" \
                --Sample_ID ${Sample_ID} \
                --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} \
                --OUTPUT_FMC_PATH ${OUTPUT_FMC_PATH} --OUTPUT_FMC_HF_PATH ${OUTPUT_FMC_HF_PATH}  --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
                --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH} \
                --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --MIN_BQ ${MIN_BQ} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} \
                --BLACKLIST ${BLACKLIST}

        fi
    done

    HOLD_J="${HOLD_J:1}"  # 맨 앞 ,를 빼줌
   

done
