#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/home/goldpm1/Meningioma/02.Align"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "signeR_01" "signeR_02"; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

REF="/home/goldpm1/reference/genome.fa"
#REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz"
    
SIGNER_PATH="/home/goldpm1/Meningioma/41.Signature/02.signeR"
UNIQUE_VCF_DIR="/home/goldpm1/Meningioma/04.mutect/05.unique"
SHARED_VCF_DIR="/home/goldpm1/Meningioma/04.mutect/03.vep"


SIGNER_INPUT_VCF_DIR=${SIGNER_PATH}"/01.vcf"
SIGNER_RESULT_DIR=${SIGNER_PATH}"/02.result"
for dir in $SIGNER_INPUT_VCF_DIR $SIGNER_RESULT_DIR ; do
    if [ ! -d ${dir} ] ; then
        mkdir -p ${dir}
    fi
done


sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬



#01. 한 곳으로 옮기기 (Sample name 2~3개 중에 1개만 뺴서 옮기기  (signeR의 특징) )
for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102
    for TISSUE in Dura Tumor; do
        if [ ! ${Sample_ID}"_"${TISSUE}  == "221102_Dura" ]; then
            echo -e ${Sample_ID}"_"${TISSUE}
            UNIQUE_VCF_PATH=${UNIQUE_VCF_DIR}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.rescue.unique.vcf"
            bcftools view -s ${Sample_ID}"_"${TISSUE} ${UNIQUE_VCF_PATH} > ${SIGNER_INPUT_VCF_DIR}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.rescue.unique.vcf"
        fi
    done
    SHARED_VCF_PATH=${SHARED_VCF_DIR}"/"${Sample_ID}"_multiple.MT2.FMC.HF.RMBLACK.vep.vcf"
    bcftools view -s ${Sample_ID}"_Tumor" ${SHARED_VCF_PATH} > ${SIGNER_INPUT_VCF_DIR}"/"${Sample_ID}"_multiple.MT2.FMC.HF.RMBLACK.vep.temp.vcf"
    # 220930_Tumor → 220930_Multiple로 바꾸기
    bcftools reheader --samples <(echo -e ${Sample_ID}_Tumor'\t'${Sample_ID}_Multiple) ${SIGNER_INPUT_VCF_DIR}"/"${Sample_ID}"_multiple.MT2.FMC.HF.RMBLACK.vep.temp.vcf" -o ${SIGNER_INPUT_VCF_DIR}"/"${Sample_ID}"_multiple.MT2.FMC.HF.RMBLACK.vep.vcf"
    rm -rf ${SIGNER_INPUT_VCF_DIR}"/"${Sample_ID}"_multiple.MT2.FMC.HF.RMBLACK.vep.temp.vcf"
done


#02. signeR 실행하기
qsub -pe smp 2 -e $logPath"/signeR_01" -o $logPath"/signeR_01" -N "signeR_01."${Sample_ID}  ${CURRENT_PATH}"/signeR_pipe_01.sh" \
    --SIGNER_INPUT_VCF_DIR ${SIGNER_INPUT_VCF_DIR} \
    --SIGNER_RESULT_DIR ${SIGNER_RESULT_DIR} 
