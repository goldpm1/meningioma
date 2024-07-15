#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/00.Targetted"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "51.cutbam"; do
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
        OUTPUT_CUTBAM_PATH=${DATA_PATH}"/09.pysam/"${hg}"/"${TISSUE}"/01.NF2_cutbam/"${Sample_ID}"_"${TISSUE}".NF2.bam"
        
        for subdir in ${OUTPUT_CUTBAM_PATH%/*} ; do
            if [ ! -d $subdir ] ; then
                mkdir -p $subdir
            fi
        done



        # 51. cutbam 수행
        qsub -pe smp 6 -e ${logPath}"/51.cutbam" -o ${logPath}"/51.cutbam" -N "cut1_"${Sample_ID}"_"${TISSUE} "06.pysam_cutbam_pipe_01.sh" \
        --REF ${REF_hg38}  --INTERVAL ${INTERVAL} --INPUT_BAM_PATH ${INPUT_BAM_PATH}  --OUTPUT_CUTBAM_PATH ${OUTPUT_CUTBAM_PATH}
        
        # 52. remove_nonref

    done
done
