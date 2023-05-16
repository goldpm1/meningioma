#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/home/goldpm1/Meningioma/02.Align"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "01.cosmic" "02.denovo"; do
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
    
SIGPROFILER_PATH="/home/goldpm1/Meningioma/41.Signature/01.SigProfiler"
UNIQUE_VCF_DIR="/home/goldpm1/Meningioma/04.mutect/05.unique"
SHARED_VCF_DIR="/home/goldpm1/Meningioma/04.mutect/03.vep"


SIGPROFILER_INPUT_VCF_DIR=${SIGPROFILER_PATH}"/01.vcf"
SIGPROFILER_RESULT_COSMIC_DIR=${SIGPROFILER_PATH}"/02.result_cosmic"
SIGPROFILER_RESULT_EXTRACT_DIR=${SIGPROFILER_PATH}"/02.result_extract"

rm -rf ${SIGPROFILER_INPUT_VCF_DIR}
if [ ! -d ${SIGPROFILER_INPUT_VCF_DIR} ] ; then
    mkdir -p ${SIGPROFILER_INPUT_VCF_DIR}
fi
if [ ! -d ${SIGPROFILER_RESULT_COSMIC_DIR} ] ; then
    mkdir -p ${SIGPROFILER_RESULT_COSMIC_DIR}
fi
if [ ! -d ${SIGPROFILER_RESULT_EXTRACT_DIR} ] ; then
    mkdir -p ${SIGPROFILER_RESULT_EXTRACT_DIR}
fi

sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬



#01. 한 곳으로 옮기기
for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102
    for TISSUE in Dura Tumor; do
        echo -e ${Sample_ID}"_"${TISSUE}
        UNIQUE_VCF_PATH=${UNIQUE_VCF_DIR}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.rescue.unique.vcf"
        SHARED_VCF_PATH=${SHARED_VCF_DIR}"/"${Sample_ID}"_multiple.MT2.FMC.HF.RMBLACK.vep.vcf"

        cp ${UNIQUE_VCF_PATH} ${SIGPROFILER_INPUT_VCF_DIR}
        cp ${SHARED_VCF_PATH} ${SIGPROFILER_INPUT_VCF_DIR}
    done
done


#02. Sigprofiler 실행하기

echo -e "qsub -pe smp 2 -e $logPath"/01.cosmic" -o $logPath"/01.cosmic" -N 'sig_01.'${Sample_ID}  ${CURRENT_PATH}"/sigprofiler_pipe_01.cosmic.sh" \
    --Sample_ID ${Sample_ID} \
    --SIGPROFILER_INPUT_VCF_DIR ${SIGPROFILER_INPUT_VCF_DIR} \
    --SIGPROFILER_RESULT_COSMIC_DIR ${SIGPROFILER_RESULT_COSMIC_DIR} \
    --SIGPROFILER_RESULT_EXTRACT_DIR ${SIGPROFILER_RESULT_EXTRACT_DIR}"

qsub -pe smp 2 -e $logPath"/01.cosmic" -o $logPath"/01.cosmic" -N 'sig_01.'${Sample_ID}  ${CURRENT_PATH}"/sigprofiler_pipe_01.cosmic.sh" \
    --Sample_ID ${Sample_ID} \
    --SIGPROFILER_INPUT_VCF_DIR ${SIGPROFILER_INPUT_VCF_DIR} \
    --SIGPROFILER_RESULT_COSMIC_DIR ${SIGPROFILER_RESULT_COSMIC_DIR} \
    --SIGPROFILER_RESULT_EXTRACT_DIR ${SIGPROFILER_RESULT_EXTRACT_DIR}
