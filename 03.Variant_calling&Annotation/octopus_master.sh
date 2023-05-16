#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/home/goldpm1/Meningioma/02.Align"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in 21.OCTcall 22.OCTpass; do
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
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102
    
    echo $Sample_ID
    
    REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
    INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
    
    for TISSUE in Tumor Dura; do
        CASE_BAM_PATH=${DATA_PATH}"/"${TISSUE}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
        CONTROL_BAM_PATH=${DATA_PATH}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
        OUTPUT_VCF_PATH=${DATA_PATH%/*}"/03.octopus/01.raw/"${Sample_ID}"_"${TISSUE}".vcf"
        
        for folder in  ${OUTPUT_VCF_PATH%/*}  ; do
            if [ ! -d $folder ] ; then
                mkdir $folder
            fi
        done
        
        #1. Octputs cancer call
        qsub -pe smp 10 -e $logPath"/21.OCTcall" -o $logPath"/21.OCTcall" -N 'OCT_01.'${Sample_ID}"_"${TISSUE} -hold_jid "doc_"${Sample_ID}"_"${TISSUE} octopus_pipe_01.call.sh \
        --Sample_ID ${Sample_ID} --CASE_BAM_PATH ${CASE_BAM_PATH} --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} --OUTPUT_VCF_PATH ${OUTPUT_VCF_PATH}  --REF ${REF}
    done
    
    #3. Multiple sample Octopus call
    CASE_BAM_PATH1=${DATA_PATH}"/Tumor/05.Final_bam/"${Sample_ID}"_Tumor.bam"
    CASE_BAM_PATH2=${DATA_PATH}"/Dura/05.Final_bam/"${Sample_ID}"_Dura.bam"
    CONTROL_BAM_PATH=${DATA_PATH}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
    OUTPUT_VCF_PATH=${DATA_PATH%/*}"/03.octopus/01.raw/"${Sample_ID}"_multiple.vcf"
    
    # qsub -pe smp 6 -e $logPath"/21.OCTcall" -o $logPath"/22.OCTcall" -N 'OCT_02.'${Sample_ID} -hold_jid "doc_"${Sample_ID}"_"${TISSUE} octopus_pipe_02.multipleOCT.sh \
    # --Sample_ID ${Sample_ID} --CASE_BAM_PATH1 ${CASE_BAM_PATH1} --CASE_BAM_PATH2 ${CASE_BAM_PATH2}  --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} \
    # --OUTPUT_VCF_PATH ${OUTPUT_VCF_PATH}  --REF ${REF}
    
done