#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"


if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "01.cosmic_vcf" "11.MatrixFormation" "12.MatrixGenerator" "14.DeNovoExtractor" "15.CosineSimilarity"; do
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

PROJECT_DIR="/data/project/Meningioma/99.Meningioma_public"
BAM_DIR=${PROJECT_DIR}"/02.bam"
MUTECT_DIR=${PROJECT_DIR}"/04.mutect/3.PASS"
SIGPROFILER_PATH=${PROJECT_DIR}"/41.Signature/01.SigProfiler"


SIGPROFILER_INPUT_VCF_DIR=${SIGPROFILER_PATH}"/01.vcf"
SIGPROFILER_RESULT_COSMIC_DIR=${SIGPROFILER_PATH}"/02.result_cosmic"
SIGPROFILER_INPUT_MATRIX_DIR=${SIGPROFILER_PATH}"/11.matrix"

if [ -d ${SIGPROFILER_INPUT_VCF_DIR} ] ; then
    rm -rf ${SIGPROFILER_INPUT_VCF_DIR}
fi
if [ -d ${SIGPROFILER_RESULT_COSMIC_DIR} ] ; then
    rm -rf ${SIGPROFILER_RESULT_COSMIC_DIR}
fi
if [ -d ${SIGPROFILER_INPUT_MATRIX_DIR} ] ; then
    rm -rf ${SIGPROFILER_INPUT_MATRIX_DIR}
fi

if [ ! -d ${SIGPROFILER_INPUT_VCF_DIR} ] ; then
    mkdir -p ${SIGPROFILER_INPUT_VCF_DIR}
fi
if [ ! -d ${SIGPROFILER_RESULT_COSMIC_DIR} ] ; then
    mkdir -p ${SIGPROFILER_RESULT_COSMIC_DIR}
fi
if [ ! -d ${SIGPROFILER_INPUT_MATRIX_DIR} ] ; then
    mkdir -p ${SIGPROFILER_INPUT_MATRIX_DIR}
fi



sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬



#01. Sample 종류별로 모으기 

for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        # MEN0045-C, MEN0045-T, ...
    VCF_PATH=${MUTECT_DIR}"/"${Sample_ID}".vcf"
    cp ${VCF_PATH} ${SIGPROFILER_INPUT_VCF_DIR}
done


######################################################## VCF : SAMPLE 단위 ################################################################3


qsub -pe smp 2 -e $logPath"/01.cosmic_vcf" -o $logPath"/01.cosmic_vcf" -N 'sig_01.cosmic_vcf'  ${CURRENT_PATH}"/sigprofiler_pipe_01.cosmic_vcf.sh" \
    --SIGPROFILER_INPUT_VCF_DIR ${SIGPROFILER_INPUT_VCF_DIR} \
    --SIGPROFILER_RESULT_COSMIC_DIR ${SIGPROFILER_RESULT_COSMIC_DIR}





######################################################## TXT : TISSUE 단위 or SAMPLE 단위 ###################################################


echo -e "\n--------------------------------------- #. SigProfiler : MatrixFormation & MatrixGenerator --------------------------"
TISSUE="all"
SIGPROFILER_INPUT_MATRIX_DIR=${SIGPROFILER_PATH}"/11.matrix/"${TISSUE}
#"BY_TISSUE" or "BY_SAMPLE"
qsub -pe smp 1 -e $logPath"/11.MatrixFormation" -o $logPath"/11.MatrixFormation" -N 'sig_11.'${TISSUE}  ${CURRENT_PATH}"/sigprofiler_pipe_11.MatrixFormation.sh" \
    --SCRIPT_DIR ${CURRENT_PATH} \
    --RUN "BY_TISSUE" \
    --TISSUE ${TISSUE} \
    --SIGPROFILER_INPUT_VCF_DIR ${SIGPROFILER_INPUT_VCF_DIR} \
    --SIGPROFILER_INPUT_MATRIX_DIR ${SIGPROFILER_INPUT_MATRIX_DIR}

    # echo -e "\n---------------------------------------- #. SigProfiler : COSMIC assignment --------------------------"
    # qsub -pe smp 1 -e $logPath"/13.assignment" -o $logPath"/13.assignment" -N 'sig_13.'${TISSUE} -hold_jid 'sig_11.'${TISSUE}  ${CURRENT_PATH}"/sigprofiler_pipe_13.COSMICassignment.sh" \
    #     --SCRIPT_DIR ${CURRENT_PATH} \
    #     --OUTPUT_SBS96 ${SIGPROFILER_INPUT_MATRIX_DIR}"/output/SBS/Meningioma.SBS96.all" \
    #     --ASSIGNMENT_DIR ${SIGPROFILER_INPUT_MATRIX_DIR}"/output/Assignment" 


    # echo -e "\n---------------------------------------- #. SigProfiler : Extractor ----------------------------------------"
    # qsub -pe smp 1 -e $logPath"/14.DeNovoExtractor" -o $logPath"/14.DeNovoExtractor" -N 'sig_14.'${TISSUE} -hold_jid 'sig_13.'${TISSUE}  ${CURRENT_PATH}"/sigprofiler_pipe_14.DeNovoExtractor.sh" \
    #     --SCRIPT_DIR ${CURRENT_PATH} \
    #     --OUTPUT_SBS96 ${SIGPROFILER_INPUT_MATRIX_DIR}"/output/SBS/Meningioma.SBS96.all" \
    #     --EXTRACTOR_DIR ${SIGPROFILER_INPUT_MATRIX_DIR}"/output/DeNovoExtractor" 

