#!/bin/bash
#$ -cwd
#$ -S /bin/bash


REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
REF_hg19="/home/goldpm1/reference/Broadhg19/Homo_sapiens_assembly19.fasta"             # 왜그런지는 모르겠지만 이 reference를 써야 돌아간다
REF_hg38="/home/goldpm1/reference/genome.fa"
hg="hg38"

INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz"

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

for sublog in "facetcnv"; do
    if [ -d $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

DATA_PATH="/data/project/Meningioma/02.Align"
FACETCNV_PATH="/data/project/Meningioma/11.cnv/5.facetcnv"

# sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
# sample_name_LIST=(${sample_name_list// / })     # array로 만듬
# for idx in ${!sample_name_LIST[@]}; do
#     Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102

for Sample_ID in 230303; do
    CONTROL_BAM_PATH=${DATA_PATH}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"

    #for TISSUE in Tumor Tumor_FFT Tumor_PCT Tumor_PFT Tumor_PP Tumor_PT Dura Ventricle Cortex; do   #Tumor Dura
    for TISSUE in Dura; do   #Tumor Dura
        CASE_BAM_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE%_*}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"

        if [ -f ${CASE_BAM_PATH} ]; then     # File이 있어야만 진행
            OUTPUT_DIR=${FACETCNV_PATH}"/"${Sample_ID}"/"${TISSUE}
            OUTPUT_PREFIX=${OUTPUT_DIR}"/"${Sample_ID}
            if [ ! -d ${OUTPUT_DIR} ]; then
                mkdir -p ${OUTPUT_DIR}
            fi
        

            module load cnv_facets
            qsub -pe smp 3 -o $logPath"/facetcnv" -e $logPath"/facetcnv" -N "fac01_"${Sample_ID}"_"${TISSUE}  -hold_jid "doc_"${Sample_ID}"_"${TISSUE} "facetcnv_pipe.sh" \
                --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} --CASE_BAM_PATH ${CASE_BAM_PATH} \
                --dbSNP ${dbSNP} --TARGETS ${INTERVAL} --OUTPUT_PREFIX ${OUTPUT_PREFIX}
        fi
    done
done


