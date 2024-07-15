#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

PROJECT_DIR="/data/project/Meningioma"

REF="/home/goldpm1/reference/genome.autosomal.fa"
hg="hg38"

GTF_FILE="/home/goldpm1/resources/gencode.v38.annotation.gtf"

for sublog in 22.quantify; do  #01.WGRS
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

FASTQ_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/00.raw"
BAM_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/02.Align"
MATRIX_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/03.matrix"
TMP_PATH=${BAM_DIR}"/temp"

FASTQ=$(find "$FASTQ_DIR" -type f -name "*.fastq.gz" )
FASTQ_LIST=(${FASTQ// / })                   # 이를 배열 (list)로 만듬

# conda activate cnvpytor

for idx in ${!FASTQ_LIST[@]}              # @ 배열의 모든 element    #! : indexing
do
        FASTQ_PATH=${FASTQ_LIST[idx]}         # idx번째의 파일명을 담아둔다
        ID=${FASTQ_LIST[idx]/"${FASTQ_DIR}/"/}  
        ID=${ID%".hifi_reads.fastq.gz"}
        # echo -e "BAM_PATH: "$BAM_PATH   # /data/project/Meningioma/01.PacBio/Isoseq_RNA/02.Align/221102_Tumor.bam
        # echo -e "ID: "$ID   # 221102_Tumor

        # 02. quantification
        OUTPUT_DIR=${MATRIX_DIR}"/02.Mando/"${ID}
        if [ ! -d ${OUTPUT_DIR} ] ; then
            mkdir -p ${OUTPUT_DIR}
        fi
        qsub -pe smp 5 -e $logPath"/22.quantify" -o $logPath"/22.quantify" -N 'PacBio_22.'${ID}  ${CURRENT_PATH}"/PacBio_Mando_pipe_02.quantify.sh" \
             --GTF_FILE ${GTF_FILE} \
             --REF ${REF} \
             --FASTQ_PATH ${FASTQ_PATH} \
             --OUTPUT_DIR ${OUTPUT_DIR}
done


