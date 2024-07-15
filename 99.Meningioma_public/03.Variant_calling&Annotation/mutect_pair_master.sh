#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

# bash /data/project/Meningioma/script/99.Meningioma_public/03.Variant_calling\&Annotation/mutect_pair_master.sh

PROJECT_DIR="/data/project/Meningioma/99.Meningioma_public/SRP261564"     # SRP050339 (single)  SRP227246 (Al-mefty) SRP261564 (Saudi, single)
FASTQ_DIR=${PROJECT_DIR}"/02.fastq"
BAM_DIR=${PROJECT_DIR}"/03.bam"
TMP_PATH=${BAM_DIR}"/temp"

hg="hg38"
if [ "${hg}" == "hg38" ]; then
    REF="/home/goldpm1/reference/genome.fa"           #REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
    dbsnp="/data/public/dbSNP/b155/GRCh38/GCF_000001405.39.re.vcf.gz"
    INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
    PON="/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
    gnomad="/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
elif [ "${hg}" == "hg19" ]; then
    REF="/home/goldpm1/reference/hg19/hg19.fa"
    dbsnp="/data/public/dbSNP/b155/GRCh37/GCF_000001405.25.re.vcf.gz"
    INTERVAL="/home/goldpm1/resources/Exon.reference.GRCh38.bed"
    PON=""
    gnomad="/data/public/GATK/gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf"
fi

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi

for sublog in 01.MTcall 02.FMC_HF_RMBLACK; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬


for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}  
    HOLD_J=""
    
    
    CASE_BAM_PATH=${BAM_DIR}"/"${hg}"/05.Final_bam/"${Sample_ID}"_MNT.bam"
    CONTROL_BAM_PATH=${BAM_DIR}"/"${hg}"/05.Final_bam/"${Sample_ID}"_MNB.bam"
    OUTPUT_VCF_GZ=${BAM_DIR%/*}"/04.mutect/01.raw/"${Sample_ID}".vcf.gz"
    OUTPUT_FMC_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}".MT2.FMC.vcf"
    OUTPUT_FMC_HF_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}".MT2.FMC.HF.vcf"
    OUTPUT_FMC_HF_RMBLACK_PATH=${BAM_DIR%/*}"/04.mutect/02.PASS/"${Sample_ID}".MT2.FMC.HF.RMBLACK.vcf"


    if [ -f ${CASE_BAM_PATH} ]; then     # File이 있어야만 진행
        SAMPLE_THRESHOLD="all"
        DP_THRESHOLD=30
        ALT_THRESHOLD=1
        REMOVE_MULTIALLELIC="True"
        PASS="True"
        REMOVE_MITOCHONDRIAL_DNA="True"
        BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
        
        for folder in  ${OUTPUT_VCF_GZ%/*}   ${OUTPUT_FMC_PATH%/*} ; do
            if [ ! -d $folder ] ; then
                mkdir $folder
            fi
        done

        #01. Mutect2 call
        qsub -pe smp 5 -e $logPath"/01.MTcall" -o $logPath"/01.MTcall" -N 'MT_01.'${Sample_ID} -hold_jid "doc_"${Sample_ID}  ${CURRENT_PATH}"/mutect_pair_pipe_01.call.sh" \
        --Sample_ID ${Sample_ID} --CASE_BAM_PATH ${CASE_BAM_PATH} --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} \
        --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ}  \
        --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH}

        #02. FMC & HF & RMBLACK & VEP
        qsub -pe smp 5 -e $logPath"/02.FMC_HF_RMBLACK" -o $logPath"/02.FMC_HF_RMBLACK" -N 'MT_02.'${Sample_ID} -hold_jid  'MT_01.'${Sample_ID}  ${CURRENT_PATH}"/mutect_pair_pipe_02.FMC_HF_RMBLACK.sh" \
            --Sample_ID ${Sample_ID} \
            --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} \
            --OUTPUT_FMC_PATH ${OUTPUT_FMC_PATH} --OUTPUT_FMC_HF_PATH ${OUTPUT_FMC_HF_PATH}  --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
            --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH} \
            --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} \
            --BLACKLIST ${BLACKLIST}

    fi
    
done


