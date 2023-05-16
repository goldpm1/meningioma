#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/00.Targetted"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in 01.fastp 02.fastqc; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

sample_name_list=$(cat ${CURRENT_PATH}"/sample_name.txt")

sample_name_LIST=(${sample_name_list// / })     # array로 만듬

# 01.QC

for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102
    echo $Sample_ID
    
    for TISSUE in Dura; do
        FASTQ_PATH_1=$DATA_PATH"/01.QC/00.raw/"$TISSUE"/"$Sample_ID"_"$TISSUE"_1.fq.gz"
        FASTQ_PATH_2=$DATA_PATH"/01.QC/00.raw/"$TISSUE"/"$Sample_ID"_"$TISSUE"_2.fq.gz"

        FASTP_PATH_1=$DATA_PATH"/01.QC/01.fastp/"$TISSUE"/"$Sample_ID"_"$TISSUE".R1.fq"
        FASTP_PATH_2=$DATA_PATH"/01.QC/01.fastp/"$TISSUE"/"$Sample_ID"_"$TISSUE".R2.fq"
        FASTQC_PATH=$DATA_PATH"/01.QC/02.fastqc/"$TISSUE
        for folder in ${FASTP_PATH_1%/*} $FASTQC_PATH; do
            if [ ! -d $folder ] ; then
                mkdir -p $folder
            fi
        done
        
        #1. FASTP
        qsub -pe smp 5 -e $logPath"/01.fastp" -o $logPath"/01.fastp" -N 'FP_'${Sample_ID}"_"${TISSUE} "01.QC_pipe_01.fastp.sh" \
        --Sample_ID ${Sample_ID}  --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2} --FASTP_PATH_1 ${FASTP_PATH_1} --FASTP_PATH_2 ${FASTP_PATH_2}

        #2. FASTQC
        qsub -pe smp 5 -e $logPath"/02.fastqc" -o $logPath"/02.fastqc" -N 'FQC_'${Sample_ID}"_"${TISSUE} -hold_jid 'FP_'${Sample_ID}"_"${TISSUE} "01.QC_pipe_02.fastqc.sh" \
        --Sample_ID ${Sample_ID}  --FASTQ_PATH_1 ${FASTP_PATH_1}".gz" --FASTQ_PATH_2 ${FASTP_PATH_2}".gz" --FASTQC_PATH ${FASTQC_PATH}

    done    
done