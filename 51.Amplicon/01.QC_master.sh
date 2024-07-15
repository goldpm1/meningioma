#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"


if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
# for sublog in 01.fastp 02.fastqc; do
#     if [ $logPath"/"$sublog ] ; then
#         rm -rf $logPath"/"$sublog
#     fi
#     if [ ! -d $logPath"/"$sublog ] ; then
#         mkdir -p $logPath"/"$sublog
#     fi
# done


################################################################# SINGLE ###################################################################
DATA_PATH="/data/project/Meningioma/51.Amplicon/01.single"

# 01.QC
for Sample_ID in 190426_NF2 220930_NF2 221026_NF2 230323_NF2 230920_NF2 230405_TRAF7; do  #
    FQ=$(find "$DATA_PATH"/01.QC/00.raw/"${Sample_ID}" -name  "*1.fq.gz"  )
    FQ_LIST=(${FQ// / })                   # 이를 배열 (list)로 만듬

    if [ ! ${#FQ_LIST[@]} -eq 0 ]; then
        echo ${Sample_ID}
    fi

    for idx in ${!FQ_LIST[@]};do              # @ 배열의 모든 element    #! : indexing
        FASTQ_PATH_1=${FQ_LIST[idx]}         # idx번째의 파일명을 담아둔다
        FASTQ_PATH_2=${FASTQ_PATH_1%"1.fq.gz"}"2.fq.gz"

        S1=${FQ_LIST[idx]##*/}
        DATE_TISSUE=${S1%_1*.fq.gz}
        DATE=${S1%_[DVC]*}

        echo -e "Single site\t"${DATE_TISSUE}

        FASTP_PATH_1=$DATA_PATH"/01.QC/01.fastp/"${Sample_ID}"/"${DATE_TISSUE}".R1.fq"
        FASTP_PATH_2=$DATA_PATH"/01.QC/01.fastp/"${Sample_ID}"/"${DATE_TISSUE}".R2.fq"
        FASTQC_PATH=$DATA_PATH"/01.QC/02.fastqc/"${Sample_ID}
        for folder in ${FASTP_PATH_1%/*} $FASTQC_PATH; do
            if [ ! -d $folder ] ; then
                mkdir -p $folder
            fi
        done
        
        #1. FASTP
        qsub -pe smp 5 -e $logPath"/01.fastp" -o $logPath"/01.fastp" -N 'FP_'${DATE_TISSUE} "01.QC_pipe_01.fastp.sh" \
        --Sample_ID ${Sample_ID}  --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2} --FASTP_PATH_1 ${FASTP_PATH_1} --FASTP_PATH_2 ${FASTP_PATH_2}

        #2. FASTQC
        qsub -pe smp 5 -e $logPath"/02.fastqc" -o $logPath"/02.fastqc" -N 'FQC_'${DATE_TISSUE} -hold_jid 'FP_'${DATE_TISSUE} "01.QC_pipe_02.fastqc.sh" \
        --Sample_ID ${Sample_ID}  --FASTQ_PATH_1 ${FASTP_PATH_1}".gz" --FASTQ_PATH_2 ${FASTP_PATH_2}".gz" --FASTQC_PATH ${FASTQC_PATH}
    
    done    
done


################################################################# MULTIPLEX  ###################################################################

DATA_PATH="/data/project/Meningioma/51.Amplicon/02.multiplex"

# 01.QC

# for Sample_ID in 221026_Dura; do
#     FASTQ_PATH_1=$(find "$DATA_PATH"/01.QC/00.raw/ -name  ${Sample_ID}*1.fq.gz  )
#     FASTQ_PATH_2=${FASTQ_PATH_1%"1.fq.gz"}"2.fq.gz"

#     S1=${FASTQ_PATH_1##*/}
#     DATE_TISSUE=${S1%_1*.fq.gz}
#     DATE=${S1%_[DVC]*}

    
#     echo -e "Multiplex\t"${DATE_TISSUE}

#     FASTP_PATH_1=$DATA_PATH"/01.QC/01.fastp/"${DATE_TISSUE}".R1.fq"
#     FASTP_PATH_2=$DATA_PATH"/01.QC/01.fastp/"${DATE_TISSUE}".R2.fq"
#     FASTQC_PATH=$DATA_PATH"/01.QC/02.fastqc/"
#     for folder in ${FASTP_PATH_1%/*} $FASTQC_PATH; do
#         if [ ! -d $folder ] ; then
#             mkdir -p $folder
#         fi
#     done
        
#     #1. FASTP
#     qsub -pe smp 5 -e $logPath"/01.fastp" -o $logPath"/01.fastp" -N 'FP_'${DATE_TISSUE} "01.QC_pipe_01.fastp.sh" \
#     --Sample_ID ${DATE_TISSUE}  --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2} --FASTP_PATH_1 ${FASTP_PATH_1} --FASTP_PATH_2 ${FASTP_PATH_2}

#     #2. FASTQC
#     qsub -pe smp 5 -e $logPath"/02.fastqc" -o $logPath"/02.fastqc" -N 'FQC_'${DATE_TISSUE} -hold_jid 'FP_'${DATE_TISSUE} "01.QC_pipe_02.fastqc.sh" \
#     --Sample_ID ${DATE_TISSUE}  --FASTQ_PATH_1 ${FASTP_PATH_1}".gz" --FASTQ_PATH_2 ${FASTP_PATH_2}".gz" --FASTQC_PATH ${FASTQC_PATH}
    
# done    