#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

PROJECT_DIR="/data/project/Meningioma"

PON="/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
REF="/home/goldpm1/reference/genome.fa"
hg="hg38"
gnomad="/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
GTF_FILE="/home/goldpm1/resources/gencode.v38.annotation.gtf"
GTF_COMPATIBLE_FILE="/home/goldpm1/tools/LIQA/gencode.v38.annotation.refgene"

GTF_FILE="/home/goldpm1/resources/gencode.v38.annotation.NF2.gtf"
GTF_COMPATIBLE_FILE="/home/goldpm1/tools/LIQA/gencode.v38.annotation.NF2.refgene"

#GTF_COMPATIBLE_FILE="/home/goldpm1/tools/LIQA/gencode.v38.annotation.CDH2.refgene"

# for sublog in 11.make_compatible 12.quantify; do  #01.WGRS
#     if [ $logPath"/"$sublog ] ; then
#         rm -rf $logPath"/"$sublog
#     fi
#     if [ ! -d $logPath"/"$sublog ] ; then
#         mkdir -p $logPath"/"$sublog
#     fi
# done

BAM_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/02.Align"
MATRIX_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/03.matrix"
TMP_PATH=${BAM_DIR}"/temp"

for subdir in ${BAM_DIR} ${MATRIX_DIR} ${TMP_PATH}; do 
    if [ ! -d $subdir ] ; then
        mkdir -p $subdir
    fi
done

BAM=$(find "$BAM_DIR" -type f -name "*.bam" | grep -v "ubam")
BAM_LIST=(${BAM// / })                   # 이를 배열 (list)로 만듬


# 01. make_compatible_GTF   한번만 하면 됨
# qsub -pe smp 5 -e $logPath"/11.make_compatible" -o $logPath"/11.make_compatible" -N 'PacBio_11.make_compatible'  ${CURRENT_PATH}"/PacBio_LIQA_pipe_01.make_compatible.sh" \
#      --GTF_FILE ${GTF_FILE} --GTF_COMPATIBLE_FILE ${GTF_COMPATIBLE_FILE}



for idx in ${!BAM_LIST[@]}              # @ 배열의 모든 element    #! : indexing
do
        BAM_PATH=${BAM_LIST[idx]}         # idx번째의 파일명을 담아둔다
        ID=${BAM_LIST[idx]/"${BAM_DIR}/"/}  
        ID=${ID%".bam"}
        # echo -e "BAM_PATH: "$BAM_PATH   # /data/project/Meningioma/01.PacBio/Isoseq_RNA/02.Align/221102_Tumor.bam
        # echo -e "ID: "$ID   # 221102_Tumor

        # 02. quantification
        OUTPUT_DIR=${MATRIX_DIR}"/01.LIQA/"${ID}
        if [ ! -d ${OUTPUT_DIR} ] ; then
            mkdir -p ${OUTPUT_DIR}
        fi
        qsub -pe smp 5 -e $logPath"/12.quantify" -o $logPath"/12.quantify" -N 'PacBio_12.'${ID}  ${CURRENT_PATH}"/PacBio_LIQA_pipe_02.quantify.sh" \
             --GTF_COMPATIBLE_FILE ${GTF_COMPATIBLE_FILE} \
             --BAM_PATH ${BAM_PATH} \
             --OUTPUT_DIR ${OUTPUT_DIR}"/"${ID}"_NF2.txt"
done
