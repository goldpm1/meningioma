#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

FASTQ_DIR="/data/project/Meningioma/71.Methylation/00.raw"
BAM_DIR="/data/project/Meningioma/71.Methylation/02.align"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in 01.align; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


       

GENOME_FOLDER="/home/goldpm1/reference/bismark_index/"

# for Sample_ID in 221026 221102 230323 230526; do 
#     for TISSUE in Tumor ; do
#         FASTQ_PATH_1=${FASTQ_DIR}"/"${Sample_ID}"_"${TISSUE}"/"${Sample_ID}"_"${TISSUE}"_1.fq.gz"
#         FASTQ_PATH_2=${FASTQ_DIR}"/"${Sample_ID}"_"${TISSUE}"/"${Sample_ID}"_"${TISSUE}"_2.fq.gz"

#         if [ -f ${FASTQ_PATH_1} ]; then     # File이 있어야만 진행
#             #echo -e ${FASTQ_PATH_1}" is exist"

#             # FASTP_PATH_1=$DATA_PATH"/01.fastp/"${TISSUE%_*}"/"$Sample_ID"_"$TISSUE".R1.fq"
#             # FASTP_PATH_2=$DATA_PATH"/01.fastp/"${TISSUE%_*}"/"$Sample_ID"_"$TISSUE".R2.fq"
#             BAM_PATH=${BAM_DIR}"/"${Sample_ID}"_"${TISSUE}
#             OUTPUT_SORT_BAM=${BAM_DIR}"/"${Sample_ID}"_"${TISSUE}"/"${Sample_ID}"_"${TISSUE}".sorted.bam"
            
#             for folder in ${BAM_PATH} ; do
#                 if [ ! -d $folder ] ; then
#                     mkdir -p $folder
#                 fi
#             done

#             qsub -pe smp 5 -e $logPath"/01.align" -o $logPath"/01.align" -N "bis_align_"${Sample_ID} "01.Bismark_align_pipe_01.sh" --GENOME_FOLDER ${GENOME_FOLDER} --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2}  --BAM_PATH ${BAM_PATH}  --OUTPUT_SORT_BAM ${OUTPUT_SORT_BAM}
#         fi
#     done
# done


#for Sample_ID in SRR13239466 SRR13239467 SRR13239468 SRR13239469 SRR13239470 SRR13239471; do 
for Sample_ID in SRR23826891 SRR23826892 SRR23826894 SRR23826896 SRR23826897; do 
        FASTQ_PATH_1=${FASTQ_DIR}"/01.public/"${Sample_ID}"/"${Sample_ID}"_1.fq.gz"
        FASTQ_PATH_2=${FASTQ_DIR}"/01.public/"${Sample_ID}"/"${Sample_ID}"_2.fq.gz"
    
        BAM_PATH=${BAM_DIR}"/01.public/"${Sample_ID}
        OUTPUT_SORT_BAM=${BAM_DIR}"/01.public/"${Sample_ID}"/"${Sample_ID}".sorted.bam"

        for folder in ${BAM_PATH} ; do
            if [ ! -d $folder ] ; then
                mkdir -p $folder
            fi
        done

        qsub -pe smp 5 -e $logPath"/01.align" -o $logPath"/01.align" -N "bis_align_"${Sample_ID} "01.Bismark_align_pipe_02.public.sh" --GENOME_FOLDER ${GENOME_FOLDER} --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2}  --BAM_PATH ${BAM_PATH}  --OUTPUT_SORT_BAM ${OUTPUT_SORT_BAM}
done
