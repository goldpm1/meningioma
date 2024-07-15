#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/00.Targetted"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "41.samtools" "42.parsing"; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done



sample_name_list=$(cat ${CURRENT_PATH}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬


for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102
    
    
    for TISSUE in Dura; do
        REF_hg38="/home/goldpm1/reference/genome.fa"
        REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
        hg="hg38"

        # INTERVAL="/home/goldpm1/resources/Exon.reference.GRCh38.bed"
        # INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
        # INTERVAL="/home/goldpm1/resources/TMB359.theragen.hg38.bed"
        INTERVAL="/home/goldpm1/resources/NF2.exome.proteincoding.bed"

        INPUT_BAM_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${TISSUE}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
        OUTPUT_MPILEUP_PATH=${DATA_PATH}"/07.mpileup/"${hg}"/"${TISSUE}"/01.mpileup/"${Sample_ID}"_"${TISSUE}".vcf.gz"
        
        for subdir in ${OUTPUT_MPILEUP_PATH%/*} ${OUTPUT_MPILEUP_CALL_PATH%/*}; do
            if [ ! -d $subdir ] ; then
                mkdir -p $subdir
            fi
        done



        # 41. samtools mpileup 수행
        qsub -pe smp 6 -e ${logPath}"/41.samtools" -o ${logPath}"/41.samtools" -N "mp1_"${Sample_ID}"_"${TISSUE} "04.mpileup_pipe_01.samtools.sh" \
        --REF ${REF_hg38}  --INTERVAL ${INTERVAL} --INPUT_BAM_PATH ${INPUT_BAM_PATH}  --OUTPUT_MPILEUP_PATH ${OUTPUT_MPILEUP_PATH}
        
        # 42. parsing 수행
        BQ="15"
        OUTPUT_MPILEUP_CALL_PATH=${DATA_PATH}"/07.mpileup/"${hg}"/"${TISSUE}"/02.parsing/"${Sample_ID}"_"${TISSUE}".mpileup_B"${BQ}".call"

    done
done
