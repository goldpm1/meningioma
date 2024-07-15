#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/61.Lowinput/01.XT_HS"
#DATA_PATH="/data/project/Meningioma/61.Lowinput/02.PTA"

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

# 01.QC

for Sample_ID in 190426; do
    for Clone_No in 3 4; do
# for Sample_ID in 230405; do
#     for Clone_No in 1 2 3 4 5 6 7 8; do
        FASTQ_PATH_1=$DATA_PATH"/01.QC/00.raw/Tumor/"$Sample_ID"_"${Clone_No}"_1.fq.gz"
        FASTQ_PATH_2=$DATA_PATH"/01.QC/00.raw/Tumor/"$Sample_ID"_"${Clone_No}"_2.fq.gz"

        echo ${FASTQ_PATH_1}

        if [ ! -f ${FASTQ_PATH_1} ] ; then
            continue  # Skip to the next iteration
        fi
        
        echo ${Sample_ID}"_"${Clone_No}

        FASTP_PATH_1=$DATA_PATH"/01.QC/01.fastp/Tumor/"$Sample_ID"_"${Clone_No}".R1.fq"
        FASTP_PATH_2=$DATA_PATH"/01.QC/01.fastp/Tumor/"$Sample_ID"_"${Clone_No}".R2.fq"
        FASTQC_PATH=$DATA_PATH"/01.QC/02.fastqc/Tumor"
        for folder in ${FASTP_PATH_1%/*} $FASTQC_PATH; do
            if [ ! -d $folder ] ; then
                mkdir -p $folder
            fi
        done
        
        #1. FASTP
        qsub -pe smp 5 -e $logPath"/01.fastp" -o $logPath"/01.fastp" -N 'FP_'${Sample_ID}"_"${Clone_No} "01.QC_pipe_01.fastp.sh" \
        --Sample_ID ${Sample_ID}  --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2} --FASTP_PATH_1 ${FASTP_PATH_1} --FASTP_PATH_2 ${FASTP_PATH_2}

        #2. FASTQC
        qsub -pe smp 5 -e $logPath"/02.fastqc" -o $logPath"/02.fastqc" -N 'FQC_'${Sample_ID}"_"${Clone_No} -hold_jid 'FP_'${Sample_ID}"_"${Clone_No} "01.QC_pipe_02.fastqc.sh" \
        --Sample_ID ${Sample_ID}  --FASTQ_PATH_1 ${FASTP_PATH_1}".gz" --FASTQ_PATH_2 ${FASTP_PATH_2}".gz" --FASTQC_PATH ${FASTQC_PATH}

    done    
done