#!/bin/bash
#$ -S /bin/bash
#$ -cwd

PROJECT_DIR="/data/project/Meningioma/61.Lowinput/01.XT_HS"
#PROJECT_DIR="/data/project/Meningioma/61.Lowinput/02.PTA"

WES_BAM_DIR="/data/project/Meningioma/02.Align"
BAM_DIR=${PROJECT_DIR}"/02.Align"
MUTECT_DIR="/data/project/Meningioma/04.mutect"
SCRIPT_DIR="/data/project/Meningioma/script/61.Lowinput"
logPath=${SCRIPT_DIR}"/log"
PON="/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
#REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
REF="/home/goldpm1/reference/genome.fa"
hg="hg38"
gnomad="/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz" 
BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
TMP_PATH=${BAM_DIR}"/temp"
INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
INTERVAL="/home/goldpm1/resources/whole.exome.bed"




for sublog in "04-1.call"; do
    if [ -d ${logPath}"/"${sublog} ]; then
        rm -rf ${logPath}"/"${sublog}
    fi
    if [ ! -d ${logPath}"/"${sublog} ]; then
        mkdir -p ${logPath}"/"${sublog}
    fi
done
for subdir in "01.call" "02.PASS"  ; do
    if [ ! -d ${MUTECT_DIR}"/"${subdir} ]; then
        mkdir -p ${MUTECT_DIR}"/"${subdir}
    fi
done



# 01. call (single)
for Sample_ID in 190426; do
#for Sample_ID in 230405; do
    BAM_DIR_LIST=""

    for Clone_No in 1 2 3 4 5 6 7 8 ; do
        CASE_BAM_PATH=${BAM_DIR}"/"${hg}"/Tumor/06.Final_bam/"${Sample_ID}"_"${Clone_No}".bam"
        BAM_DIR_LIST=${BAM_DIR_LIST}","${CASE_BAM_PATH}
        TITLE_LIST=${TITLE_LIST}","${Sample_ID}"_"${Clone_No}
        echo -e ${Sample_ID}"_"${Clone_No}
    done
      
    CONTROL_BAM_PATH=${WES_BAM_DIR}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
    BAM_DIR_LIST=${BAM_DIR_LIST}","${CONTROL_BAM_PATH}



    #01. Mutect Call
    MUTECT_OUTPUT_PATH=${MUTECT_DIR}"/01.call/"${Sample_ID}".multiple.vcf"
    MUTECT_OUTPUT_FMC_PATH=${MUTECT_DIR}"/02.PASS/"${Sample_ID}".multiple.FMC.vcf"
    MUTECT_OUTPUT_FMC_HF_PATH=${MUTECT_DIR}"/02.PASS/"${Sample_ID}".multiple.FMC.HF.vcf"
    MUTECT_OUTPUT_FMC_HF_RMBLACK_PATH=${MUTECT_DIR}"/02.PASS/"${Sample_ID}".multiple.FMC.HF.RMBLACK.vcf"

    SAMPLE_THRESHOLD="all"
    DP_THRESHOLD=30
    ALT_THRESHOLD=2
    REMOVE_MULTIALLELIC="True"
    PASS="True"
    REMOVE_MITOCHONDRIAL_DNA="True"
    BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
    
    qsub -pe smp 5 -e $logPath"/04-1.call" -o $logPath"/04-1.call" -N 'MT.01_'${Sample_ID} ${SCRIPT_DIR}"/04.mutect_multiple_pipe_01.call.sh"  \
        --SCRIPT_DIR ${SCRIPT_DIR} \
        --Sample_ID ${Sample_ID} \
        --REF ${REF} \
        --BAM_DIR_LIST ${BAM_DIR_LIST} \
        --MUTECT_OUTPUT_PATH ${MUTECT_OUTPUT_PATH} \
        --MUTECT_OUTPUT_FMC_PATH ${MUTECT_OUTPUT_PATH} \
        --PON ${PON} \
        --gnomad ${gnomad} \
        --TMP_PATH ${TMP_PATH} \
        --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} \
        --BLACKLIST ${BLACKLIST}

    # python3 ${SCRIPT_DIR}"/04.mutect_multiple_pipe_01.call.py" \
    #     --REF ${REF} \
    #     --BAM_DIR_LIST ${BAM_DIR_LIST} \
    #     --normal ${Sample_ID}"_Blood" \
    #     --panel_of_normals ${PON} \
    #     --germline_resource ${gnomad} \
    #     --O ${MUTECT_OUTPUT_VCF} \
    #     --temp_dir ${TMP_PATH}

done