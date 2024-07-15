#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/00.Targetted"
WES_BAM_DIR="/data/project/Meningioma/02.Align"
BAM_DIR=${DATA_PATH}"/02.Align"
PON="/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
#REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
REF="/home/goldpm1/reference/genome.fa"
hg="hg38"
gnomad="/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
TMP_PATH=${DATA_PATH}"/temp"
INTERVAL="/home/goldpm1/resources/TMB359.theragen.hg38.bed"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi

for sublog in 11.multipleMT 12.vep; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


# sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
# sample_name_LIST=(${sample_name_list// / })     # array로 만듬
# for idx in ${!sample_name_LIST[@]}; do
#     Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102, 230127, 230323, 230419
# done


######################################  MULTIPLE CALL ########################################
for Sample_ID in 230526; do
    CASE_BAM_PATH_DURA1=${WES_BAM_DIR}"/"${hg}"/Dura/05.Final_bam/"${Sample_ID}"_Dura.bam"
    CASE_BAM_PATH_DURA2=${BAM_DIR}"/"${hg}"/Dura/06.Final_bam/"${Sample_ID}"_2_Dura.bam"
    CASE_BAM_PATH_DURA3=${BAM_DIR}"/"${hg}"/Dura/06.Final_bam/"${Sample_ID}"_3_Dura.bam"
    
    CONTROL_BAM_PATH=${WES_BAM_DIR}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"

    
    SAMPLE_THRESHOLD="all"
    DP_THRESHOLD=30
    ALT_THRESHOLD=2
    REMOVE_MULTIALLELIC="True"
    PASS="True"
    REMOVE_MITOCHONDRIAL_DNA="True"
    BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
    

    ####################################### 3개만 있는 경우 (Tumor, Dura, Cortex) #######################################

    for NUM in "12" "13" "23"; do
        if [[ ${NUM} == "12" ]]; then
            CASE_BAM_PATH1=${CASE_BAM_PATH_DURA1}
            CASE_BAM_PATH2=${CASE_BAM_PATH_DURA2}
        fi
        if [[ ${NUM} == "13" ]]; then
            CASE_BAM_PATH1=${CASE_BAM_PATH_DURA1}
            CASE_BAM_PATH2=${CASE_BAM_PATH_DURA3}
        fi
        if [[ ${NUM} == "23" ]]; then
            CASE_BAM_PATH1=${CASE_BAM_PATH_DURA2}
            CASE_BAM_PATH2=${CASE_BAM_PATH_DURA3}
        fi

        OUTPUT_VCF_GZ=${DATA_PATH}"/11.mutect_multiple/01.raw/"${Sample_ID}"_"${NUM}".vcf.gz"
        OUTPUT_FMC_PATH=${DATA_PATH}"/11.mutect_multiple/02.PASS/"${Sample_ID}"_"${NUM}".MT2.FMC.vcf"
        OUTPUT_FMC_HF_PATH=${DATA_PATH}"/11.mutect_multiple/02.PASS/"${Sample_ID}"_"${NUM}".MT2.FMC.HF.vcf"
        OUTPUT_FMC_HF_RMBLACK_PATH=${DATA_PATH}"/11.mutect_multiple/02.PASS/"${Sample_ID}"_"${NUM}".MT2.FMC.HF.RMBLACK.vcf"

        for subdir in ${OUTPUT_VCF_GZ} ${OUTPUT_FMC_PATH} ${OUTPUT_FMC_HF_RMBLACK_PATH}; do
            if [ ! -d ${subdir%/*} ] ; then
                mkdir -p ${subdir%/*}
            fi
        done

        #11. Multiple sample Mutect2 call (FMC, HF, RMBLACK 포함)
        qsub -pe smp 8 -e $logPath"/11.multipleMT" -o $logPath"/11.multipleMT" -N "MT_11."${Sample_ID}"_"${NUM}   ${CURRENT_PATH}"/03.mutect_pipe_11.multipleMT.sh" \
        --Sample_ID ${Sample_ID} \
        --CASE_BAM_PATH1 ${CASE_BAM_PATH1} --CASE_BAM_PATH2 ${CASE_BAM_PATH2}  --NUM ${NUM}  \
        --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} \
        --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} --OUTPUT_FMC_PATH ${OUTPUT_FMC_PATH}  --OUTPUT_FMC_HF_PATH ${OUTPUT_FMC_HF_PATH} --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
        --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH} \
        --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} \
        --BLACKLIST ${BLACKLIST}

        #12. VEP annotation  +  Nearest gene annotation
        INPUT_VCF=${OUTPUT_FMC_HF_RMBLACK_PATH}
        VEP_FMC_HF_RMBLACK_PATH=${DATA_PATH}"/11.mutect_multiple/03.vep/"${Sample_ID}"_"${NUM}".MT2.FMC.HF.RMBLACK.vep.vcf"
        for subdir in ${VEP_FMC_HF_RMBLACK_PATH}; do
            if [ ! -d ${subdir%/*} ] ; then
                mkdir -p ${subdir%/*}
            fi
        done
        qsub -pe smp 6 -e $logPath"/12.vep" -o $logPath"/12.vep" -N "MT_12."${Sample_ID}"_"${NUM} -hold_jid "MT_11."${Sample_ID}"_"${NUM}  ${CURRENT_PATH}"/03.mutect_pipe_20.vep.sh" \
            --REF ${REF} --INPUT_VCF ${INPUT_VCF} --OUTPUT_VCF ${VEP_FMC_HF_RMBLACK_PATH}
    done


done