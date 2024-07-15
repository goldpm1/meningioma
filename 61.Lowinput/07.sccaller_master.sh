#!/bin/bash
#$ -S /bin/bash
#$ -cwd

PROJECT_DIR="/data/project/Meningioma/61.Lowinput/01.XT_HS"
BAM_DIR=${PROJECT_DIR}"/02.Align"
SCCALLER_DIR=${PROJECT_DIR}"/07.sccaller"
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




for sublog in "07-1.bulk_found" "07-2.bulk_notfound"; do
    # if [ -d ${logPath}"/"${sublog} ]; then
    #     rm -rf ${logPath}"/"${sublog}
    # fi
    if [ ! -d ${logPath}"/"${sublog} ]; then
        mkdir -p ${logPath}"/"${sublog}
    fi
done
for subdir in "07-1.bulk_found" "07-2.bulk_notfound"; do
    if [ ! -d ${SCCALLER_DIR}"/"${subdir} ]; then
        mkdir -p ${SCCALLER_DIR}"/"${subdir}
    fi
done



# 01. call (single)
for Sample_ID in 190426; do
    HC_VCF_LIST=""
    BAM_DIR_LIST=""
    TITLE_LIST=""
    hold_j=""

    WES_TUMOR_BAM_PATH="/data/project/Meningioma/02.Align/"${hg}"/Tumor/05.Final_bam/"${Sample_ID}"_Tumor.bam"
    WES_TUMOR_BED_PATH="/data/project/Meningioma/04.mutect/03.Tumor_interval/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.bed"
    if [ ${Sample_ID} == "190426" ]; then
        WES_TUMOR_BAM_PATH="/data/project/Meningioma/02.Align/"${hg}"/Tumor/05.Final_bam/"${Sample_ID}"_Tumor_PT.bam"
    fi
    CONTROL_BAM_PATH=${BAM_DIR}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"

    # bamsnap 찍기 위해 BAM_DIR_LIST 정리해줌
    for Clone_No in 1 2; do
        CASE_BAM_PATH=${BAM_DIR}"/"${hg}"/Tumor/06.Final_bam/"${Sample_ID}"_"${Clone_No}".bam"
        BAM_DIR_LIST=${BAM_DIR_LIST}","${CASE_BAM_PATH}
        TITLE_LIST=${TITLE_LIST}","${Sample_ID}"_"${Clone_No}
    done
    BAM_DIR_LIST=${BAM_DIR_LIST}","${WES_TUMOR_BAM_PATH}","${CONTROL_BAM_PATH}
    TITLE_LIST=${TITLE_LIST}","${Sample_ID}"_Tumor,"${Sample_ID}"_Blood"


    for Clone_No in 1 2; do
        CASE_BAM_PATH=${BAM_DIR}"/"${hg}"/Tumor/06.Final_bam/"${Sample_ID}"_"${Clone_No}".bam"
        echo -e ${Sample_ID}"_"${Clone_No}
            
        # 01. bulk에서 찾은 것을 sc bam에서 검증 (간단해 보이지만 1시간 걸림)
        # qsub -pe smp 2 -e $logPath"/07-1.bulk_found" -o $logPath"/07-1.bulk_found"  -N 'SCCALL1_'${Sample_ID}"_"${Clone_No} ${SCRIPT_DIR}"/07.sccaller_pipe_01.bulk_found.sh" \
        #     --CASE_BAM_PATH ${CASE_BAM_PATH} \
        #     --WES_TUMOR_BED_PATH ${WES_TUMOR_BED_PATH} \
        #     --OUTPUT_BAMSNAP_DIR ${OUTPUT_BAMSNAP_DIR} \
        #     --SAMPLE_ID ${Sample_ID} \
        #     --BAM_DIR_LIST ${BAM_DIR_LIST} \
        #     --TITLE_LIST ${TITLE_LIST} \
        #     --REF ${REF} \
        #     --OUTPUT_SCCALLER ${SCCALLER_DIR}"/07-1.bulk_found/"${Sample_ID}"_"${Clone_No}".sccaller.vcf"

        #02. bulk에서 못 찾은것을 새로 찾아주기
        qsub -pe smp 2 -e $logPath"/07-2.bulk_notfound" -o $logPath"/07-2.bulk_notfound"  -N 'SCCALL2_'${Sample_ID}"_"${Clone_No} ${SCRIPT_DIR}"/07.sccaller_pipe_02.bulk_notfound.sh" \
            --CASE_BAM_PATH ${CASE_BAM_PATH} \
            --WES_TUMOR_BAM_PATH ${WES_TUMOR_BAM_PATH} \
            --WES_TUMOR_BED_PATH ${WES_TUMOR_BED_PATH} \
            --OUTPUT_BAMSNAP_DIR ${OUTPUT_BAMSNAP_DIR} \
            --SAMPLE_ID ${Sample_ID} \
            --BAM_DIR_LIST ${BAM_DIR_LIST} \
            --TITLE_LIST ${TITLE_LIST} \
            --REF ${REF} \
            --OUTPUT_SCCALLER ${SCCALLER_DIR}"/07-2.bulk_notfound/"${Sample_ID}"_"${Clone_No}".sccaller.vcf"
    done

done


