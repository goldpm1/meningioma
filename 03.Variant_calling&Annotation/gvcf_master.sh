#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/02.Align"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "21.gvcf_call" "22.gvcf_remove_nonref"; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


GVCF_PATH="/data/project/Meningioma/05.gvcf"

sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬


for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102

    #REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
    REF="/home/goldpm1/reference/genome.fa"
    hg="hg38"
    INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
    
    for TISSUE in Tumor Dura Ventricle Cortex ; do   #Tumor Dura
        CASE_BAM_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
        CONTROL_BAM_PATH=${DATA_PATH}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"

        if [ -f ${CASE_BAM_PATH} ]; then     # File이 있어야만 진행
            echo $Sample_ID"_"$TISSUE

            #01. GVCF call
            OUTPUT_GVCF=${GVCF_PATH}"/01.call/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".g.vcf"
            if [ ! -d ${OUTPUT_GVCF%/*} ] ; then
                mkdir -p ${OUTPUT_GVCF%/*}
            fi
            qsub -pe smp 5 -e $logPath"/21.gvcf_call" -o $logPath"/21.gvcf_call" -N 'gvcf_21.'${Sample_ID}"_"${TISSUE} -hold_jid "doc_"${Sample_ID}"_"${TISSUE}  gvcf_pipe_01.call.sh \
                --CASE_BAM_PATH ${CASE_BAM_PATH} --INTERVAL ${INTERVAL} --REF ${REF} --OUTPUT_GVCF ${OUTPUT_GVCF}

            #02. remove nonref
            INPUT_GVCF=${OUTPUT_GVCF}
            OUTPUT_GVCF=${GVCF_PATH}"/02.remove_nonref/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".g.vcf"
            if [ ! -d ${OUTPUT_GVCF%/*} ] ; then
                mkdir -p ${OUTPUT_GVCF%/*}
            fi
            qsub -pe smp 1 -e $logPath"/22.gvcf_remove_nonref" -o $logPath"/22.gvcf_remove_nonref" -N 'gvcf_22.'${Sample_ID}"_"${TISSUE} -hold_jid 'gvcf_21.'${Sample_ID}"_"${TISSUE}  gvcf_pipe_02.remove_nonref.sh \
                --INPUT_GVCF ${INPUT_GVCF} --OUTPUT_GVCF ${OUTPUT_GVCF}

        fi 
    done
done