#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

#PATH#
TMP_PATH="/data/project/Meningioma/03.mpileup/hg38/Dura/temp"
#REF="/data/resource/reference/human/NCBI/GRCh38_GATK/BWAIndex/genome.fa"
REF="/home/goldpm1/reference/genome.fa"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi

for sublog in 41.mpileup 42.parsing ; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


INTERVAL="/data/project/Meningioma/03.mpileup/hg38/interval.bed"
INTERVAL="/data/project/Meningioma/03.mpileup/hg38/NF2.bed"

for date in 230526; do
#for date in 230323_2 230405_2; do
    for TISSUE in Dura Tumor; do
        echo ${date}"_"${TISSUE}

        INPUT_PATH="/data/project/Meningioma/02.Align/hg38/"${TISSUE}"/05.Final_bam/"${date}"_"${TISSUE}".bam"
        OUTPUT_PATH="/data/project/Meningioma/03.mpileup/hg38/"${TISSUE}"/01.mpileup/"${date}"_WES.pileup"

        qsub -pe smp 2 -e $logPath"/41.mpileup" -o $logPath"/41.mpileup" -N 'mpileup_'${date}"_"${TISSUE} samtoolsmpileup_231202_pipe.sh \
            --REF ${REF}  --TMP_PATH ${TMP_PATH} --INPUT_PATH ${INPUT_PATH} --INTERVAL ${INTERVAL} --OUTPUT_PATH ${OUTPUT_PATH}
        #echo "1.WES : /opt/Yonsei/samtools/1.7/samtools mpileup -f ${REF} -l ${INTERVAL} -o ${OUTPUT_PATH} ${INPUT_PATH}"
    done
    echo -e "\n"
done
