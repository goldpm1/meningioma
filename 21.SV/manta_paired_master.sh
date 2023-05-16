#!/bin/bash
#$ -cwd
#$ -S /bin/bash


############### Germline Joint Diploid ##################


REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
REF_hg19="/home/goldpm1/reference/Broadhg19/Homo_sapiens_assembly19.fasta"             # 왜그런지는 모르겠지만 이 reference를 써야 돌아간다
REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"

INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed.gz"
dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz"

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

for sublog in "02.manta_paired"; do
    if [ -d $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

DATA_PATH="/home/goldpm1/Meningioma/02.Align"
MANTA_PATH="/home/goldpm1/Meningioma/21.SV/02.manta_paired"

for subpath in "01.raw" "02.PASS"; do
    if [ ! -d $MANTA_PATH"/"$subpath ] ; then
        mkdir -p $MANTA_PATH"/"$subpath
    fi
done

sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬


for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102

    CONTROL_BAM_PATH=${DATA_PATH}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
    TISSUE="Tumor"
    TUMOR_BAM_PATH=${DATA_PATH}"/"${TISSUE}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
    TISSUE="Dura"
    DURA_BAM_PATH=${DATA_PATH}"/"${TISSUE}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"

    OUTPUT_PATH_01=${MANTA_PATH}"/01.raw/"${Sample_ID}
    OUTPUT_PATH_02=${MANTA_PATH}"/02.PASS/"${Sample_ID}".Manta.PASS.vcf"
    if [ -d ${OUTPUT_PATH_01} ] ; then
        rm -rf ${OUTPUT_PATH_01}
    fi
    if [ ! -d ${OUTPUT_PATH_01} ] ; then
        mkdir -p ${OUTPUT_PATH_01}
    fi
    
    qsub -pe smp 3 -o $logPath"/02.manta_paired" -e $logPath"/02.manta_paired" -N "manp01_"${Sample_ID}  "manta_paired_pipe_01.sh" \
        --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} --TUMOR_BAM_PATH ${TUMOR_BAM_PATH}  --DURA_BAM_PATH ${DURA_BAM_PATH} \
        --REF ${REF_hg38} --CALLREGIONS ${INTERVAL} --OUTPUT_PATH_01 ${OUTPUT_PATH_01} --OUTPUT_PATH_02 ${OUTPUT_PATH_02}
    
done


