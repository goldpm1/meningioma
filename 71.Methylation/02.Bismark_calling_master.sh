#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

FASTQ_DIR="/data/project/Meningioma/71.Methylation/00.raw"
BAM_DIR="/data/project/Meningioma/71.Methylation/02.align"
CALL_DIR_PARENT="/data/project/Meningioma/71.Methylation/03.call"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in 02.call; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


       

GENOME_FOLDER="/home/goldpm1/reference/bismark_index/"

# for Sample_ID in 221102; do 
#     for TISSUE in Tumor ; do
#         BAM_DIR2="${BAM_DIR}/${Sample_ID}_${TISSUE}"
#         BAM_PATH=$(find "${BAM_DIR2}" -name "*.bam" | head -n 1)

#         CALL_DIR=${CALL_DIR_PARENT}"/"${Sample_ID}"_"${TISSUE}

#         for folder in ${CALL_DIR} ; do
#             if [ ! -d $folder ] ; then
#                 mkdir -p $folder
#             fi
#         done

#         qsub -pe smp 5 -e $logPath"/02.call" -o $logPath"/02.call" -N "bis_call01_"${Sample_ID} "02.Bismark_calling_pipe_01.sh" \
#             --GENOME_FOLDER ${GENOME_FOLDER}  --BAM_PATH ${BAM_PATH}   --CALL_DIR ${CALL_DIR}

#         qsub -pe smp 1 -e $logPath"/02.call" -o $logPath"/02.call" -N "bis_call02_"${Sample_ID} -hold_jid "bis_call01_"${Sample_ID} "02.Bismark_calling_pipe_02.sh" \
#             --CALL_DIR ${CALL_DIR} \
#             --DSS_PATH ${CALL_DIR}"/"${Sample_ID}"_"${TISSUE}".DSS.txt"
#     done
# done



#for Sample_ID in SRR23826890 SRR23826893 SRR23826895 SRR23826898 SRR23826899; do 
for Sample_ID in SRR23826891 SRR23826892 SRR23826894 SRR23826896 SRR23826897; do 
        BAM_DIR2="${BAM_DIR}/01.public/${Sample_ID}"
        BAM_PATH=$(find "${BAM_DIR2}" -name "*.bam" ! -name "*sorted.bam"  | head -n 1)

        CALL_DIR=${CALL_DIR_PARENT}"/01.public/"${Sample_ID}

        for folder in ${CALL_DIR} ; do
            if [ ! -d $folder ] ; then
                mkdir -p $folder
            fi
        done

        qsub -pe smp 5 -e $logPath"/02.call" -o $logPath"/02.call" -hold_jid  "bis_align_"${Sample_ID} -N "bis_call01_"${Sample_ID} "02.Bismark_calling_pipe_01.sh" \
            --GENOME_FOLDER ${GENOME_FOLDER}  --BAM_PATH ${BAM_PATH}   --CALL_DIR ${CALL_DIR}

        qsub -pe smp 1 -e $logPath"/02.call" -o $logPath"/02.call"  -hold_jid "bis_call01_"${Sample_ID} -N "bis_call02_"${Sample_ID}  "02.Bismark_calling_pipe_02.sh" \
            --CALL_DIR ${CALL_DIR} \
            --DSS_PATH ${CALL_DIR}"/"${Sample_ID}".DSS.txt"
done
