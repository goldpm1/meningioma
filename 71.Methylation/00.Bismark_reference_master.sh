#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/71.Methylation"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in 00.reference; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


       
#1. Bismark index formation
GENOME_FOLDER="/home/goldpm1/reference/bismark_index"

qsub -pe smp 5 -e $logPath"/00.reference" -o $logPath"/00.reference" -N "bis_reference" "00.Bismark_reference_pipe_01.sh" --GENOME_FOLDER ${GENOME_FOLDER}

