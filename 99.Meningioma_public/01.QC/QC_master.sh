
#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

FASTQ_DIR="/data/project/Meningioma/99.Meningioma_public/SRP227246/02.fastq"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in 01.gz 02.fastp 03.fastqc; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬


for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102, 230127, 230323, 230419
    echo $Sample_ID
    
    FASTQ_RAW_PATH_1=${FASTQ_DIR}"/01.gz/"${Sample_ID}"_1.fastq"
    FASTQ_RAW_PATH_2=${FASTQ_DIR}"/01.gz/"${Sample_ID}"_2.fastq"
    FASTQ_PATH_1=${FASTQ_DIR}"/01.gz/"${Sample_ID}"_1.fq.gz"
    FASTQ_PATH_2=${FASTQ_DIR}"/01.gz"${Sample_ID}"_2.fq.gz"
    FASTP_PATH_1=${FASTQ_DIR}"/02.fastp/"${Sample_ID}".R1.fq"
    FASTP_PATH_2=${FASTQ_DIR}"/02.fastp/"${Sample_ID}".R2.fq"
    FASTQC_DIR=${FASTQ_DIR}"/03.fastqc/"$TISSUE

    for folder in ${FASTP_PATH_1%/*} $FASTQC_DIR; do
        if [ ! -d $folder ] ; then
            mkdir -p $folder
        fi
    done

    #1. gz
    qsub -pe smp 1 -e $logPath"/01.gz" -o $logPath"/01.gz" -N 'GZ_'${Sample_ID} QC_pipe_01.gz.sh \
    --Sample_ID ${Sample_ID}  --FASTQ_RAW_PATH_1 ${FASTQ_RAW_PATH_1} --FASTQ_RAW_PATH_2 ${FASTQ_RAW_PATH_2} --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2}

    #2. FASTP
    qsub -pe smp 3 -e $logPath"/02.fastp" -o $logPath"/02.fastp" -N 'FP_'${Sample_ID} -hold_jid 'GZ_'${Sample_ID} QC_pipe_02.fastp.sh \
    --Sample_ID ${Sample_ID}  --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2} --FASTP_PATH_1 ${FASTP_PATH_1} --FASTP_PATH_2 ${FASTP_PATH_2}

    #3. FASTQC
    qsub -pe smp 3 -e $logPath"/03.fastqc" -o $logPath"/03.fastqc" -N 'FQC_'${Sample_ID} -hold_jid 'FP_'${Sample_ID} QC_pipe_03.fastqc.sh \
    --Sample_ID ${Sample_ID}  --FASTQ_PATH_1 ${FASTP_PATH_1}".gz" --FASTQ_PATH_2 ${FASTP_PATH_2}".gz" --FASTQC_DIR ${FASTQC_DIR}
    
done