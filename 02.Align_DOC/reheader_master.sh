#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/02.Align"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in 00.reheader; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

REF_hg38="/home/goldpm1/reference/genome.fa"
#REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"

hg="hg38"
REF=${REF_hg38}





for Sample_ID in 230802_2 240403_2 240412_2 ; do    # 230802 240403 240412
    for TISSUE in Dura; do  

        FINAL_BAM_PATH_PREVIOUS=${DATA_PATH}"/"${hg}"/"${TISSUE%_*}"/05.Final_bam/"${Sample_ID%_*}"_"${TISSUE}".bam"
        FINAL_BAM_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE%_*}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"

        qsub -pe smp 1 -e ${logPath}"/00.reheader" -o ${logPath}"/00.reheader" -N "reheader_"${Sample_ID}"_"${TISSUE} reheader_pipe_01.sh \
            --Sample_ID ${Sample_ID} --TISSUE ${TISSUE}  --FINAL_BAM_PATH ${FINAL_BAM_PATH}  --FINAL_BAM_PATH_PREVIOUS ${FINAL_BAM_PATH_PREVIOUS}

    done
done