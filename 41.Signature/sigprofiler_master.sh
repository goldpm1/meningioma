#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/home/goldpm1/Meningioma/02.Align"

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

PROJECT_DIR="/data/project//Meningioma"
SIGPROFILER_PATH=${PROJECT_DIR}"/41.Signature/01.SigProfiler"
SHARED_VCF_DIR=${PROJECT_DIR}"/04.mutect/05.Shared_variant"
UNIQUE_VCF_DIR=${PROJECT_DIR}"/04.mutect/06.Unique"


SIGPROFILER_INPUT_VCF_DIR=${SIGPROFILER_PATH}"/01.vcf"
SIGPROFILER_RESULT_COSMIC_DIR=${SIGPROFILER_PATH}"/02.result_cosmic"
SIGPROFILER_RESULT_EXTRACT_DIR=${SIGPROFILER_PATH}"/02.result_extract"
SIGPROFILER_INPUT_MATRIX_DIR=${SIGPROFILER_PATH}"/11.matrix"

rm -rf ${SIGPROFILER_INPUT_VCF_DIR}

for TISSUE in Shared Tumor Dura; do
    if [ -d ${SIGPROFILER_INPUT_VCF_DIR}"/"${TISSUE} ] ; then
        rm -rf ${SIGPROFILER_INPUT_VCF_DIR}"/"${TISSUE}
    fi
    if [ -d ${SIGPROFILER_RESULT_COSMIC_DIR}"/"${TISSUE} ] ; then
        rm -rf ${SIGPROFILER_RESULT_COSMIC_DIR}"/"${TISSUE}
    fi
    if [ -d ${SIGPROFILER_RESULT_EXTRACT_DIR}"/"${TISSUE} ] ; then
        rm -rf ${SIGPROFILER_RESULT_EXTRACT_DIR}"/"${TISSUE}
    fi
    if [ -d ${SIGPROFILER_INPUT_MATRIX_DIR}"/"${TISSUE} ] ; then
        rm -rf ${SIGPROFILER_INPUT_MATRIX_DIR}"/"${TISSUE}
    fi
    if [ ! -d ${SIGPROFILER_INPUT_VCF_DIR}"/"${TISSUE} ] ; then
        mkdir -p ${SIGPROFILER_INPUT_VCF_DIR}"/"${TISSUE}
    fi
    if [ ! -d ${SIGPROFILER_RESULT_COSMIC_DIR}"/"${TISSUE} ] ; then
        mkdir -p ${SIGPROFILER_RESULT_COSMIC_DIR}"/"${TISSUE}
    fi
    if [ ! -d ${SIGPROFILER_RESULT_EXTRACT_DIR}"/"${TISSUE} ] ; then
        mkdir -p ${SIGPROFILER_RESULT_EXTRACT_DIR}"/"${TISSUE}
    fi
    if [ ! -d ${SIGPROFILER_INPUT_MATRIX_DIR}"/"${TISSUE} ] ; then
        mkdir -p ${SIGPROFILER_INPUT_MATRIX_DIR}"/"${TISSUE}
    fi
done


# sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
# sample_name_LIST=(${sample_name_list// / })     # array로 만듬
# for idx in ${!sample_name_LIST[@]}; do
#     Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102



#01. Sample 종류별로 모으기 
for Sample_ID in 220930 221026 221102 221202 230127 230323_2 230405_2 230419 230526 230822 230920; do
    SHARED_VCF_PATH=${SHARED_VCF_DIR}"/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.shared_variant.vcf"
    UNIQUE_VCF_PATH=${UNIQUE_VCF_DIR}"/"${Sample_ID}"_Dura.MT2.FMC.HF.RMBLACK.tumor_unique.vcf"
    cp ${SHARED_VCF_PATH} ${SIGPROFILER_INPUT_VCF_DIR}"/Shared/"${Sample_ID}"_Shared.MT2.FMC.HF.RMBLACK.shared_variant.vcf"
    cp ${UNIQUE_VCF_PATH} ${SIGPROFILER_INPUT_VCF_DIR}"/Tumor/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.tumor_unique.vcf"
    for TISSUE in Dura Ventricle Cortex; do
    UNIQUE_VCF_PATH=${UNIQUE_VCF_DIR}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.other_unique.vcf"
        if [ -f ${UNIQUE_VCF_PATH} ]; then     # File이 있어야만 진행
            cp ${UNIQUE_VCF_PATH} ${SIGPROFILER_INPUT_VCF_DIR}"/Dura"
        fi
    done
done


######################################################## 01. VCF : SAMPLE 단위 ################################################################3

for TISSUE in Shared Tumor Dura; do
    qsub -pe smp 2 -e $logPath"/01.cosmic_vcf" -o $logPath"/01.cosmic_vcf" -N 'sig_01.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/sigprofiler_pipe_01.cosmic_vcf.sh" \
        --SIGPROFILER_INPUT_VCF_DIR ${SIGPROFILER_INPUT_VCF_DIR}"/"${TISSUE} \
        --SIGPROFILER_RESULT_COSMIC_DIR ${SIGPROFILER_RESULT_COSMIC_DIR}"/"${TISSUE} \
        --SIGPROFILER_RESULT_EXTRACT_DIR ${SIGPROFILER_RESULT_EXTRACT_DIR}"/"${TISSUE}
done




# ######################################################## 11~14. TXT : TISSUE 단위 or SAMPLE 단위 ###################################################

for TISSUE in Shared Tumor Dura; do
    echo -e "\n--------------------------------------- #. SigProfiler : MatrixFormation & MatrixGenerator --------------------------"
    SIGPROFILER_INPUT_MATRIX_DIR=${SIGPROFILER_PATH}"/11.matrix/"${TISSUE}
    #"BY_TISSUE" or "BY_SAMPLE"
    qsub -pe smp 1 -e $logPath"/11.MatrixFormation" -o $logPath"/11.MatrixFormation" -N 'sig_11.'${TISSUE}  ${CURRENT_PATH}"/sigprofiler_pipe_11.MatrixFormation.sh" \
        --SCRIPT_DIR ${CURRENT_PATH} \
        --RUN "BY_TISSUE" \
        --TISSUE ${TISSUE} \
        --SIGPROFILER_INPUT_VCF_DIR ${SIGPROFILER_INPUT_VCF_DIR}"/"${TISSUE} \
        --SIGPROFILER_INPUT_MATRIX_DIR ${SIGPROFILER_INPUT_MATRIX_DIR}

    echo -e "\n---------------------------------------- #. SigProfiler : COSMIC assignment --------------------------"
    qsub -pe smp 1 -e $logPath"/13.assignment" -o $logPath"/13.assignment" -N 'sig_13.'${TISSUE} -hold_jid 'sig_11.'${TISSUE}  ${CURRENT_PATH}"/sigprofiler_pipe_13.COSMICassignment.sh" \
        --SCRIPT_DIR ${CURRENT_PATH} \
        --OUTPUT_SBS96 ${SIGPROFILER_INPUT_MATRIX_DIR}"/output/SBS/Meningioma.SBS96.all" \
        --ASSIGNMENT_DIR ${SIGPROFILER_INPUT_MATRIX_DIR}"/output/Assignment" 


    echo -e "\n---------------------------------------- #. SigProfiler : Extractor ----------------------------------------"
    qsub -pe smp 1 -e $logPath"/14.DeNovoExtractor" -o $logPath"/14.DeNovoExtractor" -N 'sig_14.'${TISSUE} -hold_jid 'sig_13.'${TISSUE}  ${CURRENT_PATH}"/sigprofiler_pipe_14.DeNovoExtractor.sh" \
        --SCRIPT_DIR ${CURRENT_PATH} \
        --OUTPUT_SBS96 ${SIGPROFILER_INPUT_MATRIX_DIR}"/output/SBS/Meningioma.SBS96.all" \
        --EXTRACTOR_DIR ${SIGPROFILER_INPUT_MATRIX_DIR}"/output/DeNovoExtractor" 



# # Cosine similiarity
# echo -e "\n---------------------------------------- #5. Cosine Similiarity ----------------------------------------"
# echo -e python3 ${SCRIPT_DIR}"/3.BioData_pipe1-5.cosinesim.py" \
#  --SIGPROFILER_PATH ${OUTPUT_DIR}"/output//DeNovo/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt" \
#  --OUTPUT_PATH ${OUTPUT_DIR}"/output/DeNovo/Cosine_sim.txt" 

# python3 ${SCRIPT_DIR}"/3.BioData_pipe1-5.cosinesim.py" \
#  --SIGPROFILER_PATH ${OUTPUT_DIR}"/output//DeNovo/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt" \
#  --OUTPUT_PATH ${OUTPUT_DIR}"/output/DeNovo/Cosine_sim.txt" 


done