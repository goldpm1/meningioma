#!/bin/bash
#$ -cwd
#$ -S /bin/bash


REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
REF_hg19="/home/goldpm1/reference/Broadhg19/Homo_sapiens_assembly19.fasta"             # 왜그런지는 모르겠지만 이 reference를 써야 돌아간다
REF_hg38="/home/goldpm1/reference/genome.fa"
hg="hg38"

INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz"

DIR="/home/goldpm1/Meningioma"

CURRENT_PATH=`pwd -P`
CURRENT_PATH='/home/goldpm1/Meningioma/script/11.cnv'
logPath=$CURRENT_PATH"/log"

DATA_PATH="/home/goldpm1/Meningioma/02.Align"
HATCHET_PATH="/home/goldpm1/Meningioma/11.cnv/4.hatchet_v1.0"

###### conda activate cnvpytor 필요


for sublog in "hatchet_v1.0"; do
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
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102

    CONTROL_BAM_PATH=${DATA_PATH}"/"${hg38}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
    CASE_BAM_PATH=""
    SAMPLES=""
    ALLNAMES=${Sample_ID}"_Blood"

    for TISSUE in Tumor Dura ; do   #Tumor Dura
        SAMPLES=${SAMPLES}","${Sample_ID}"_"${TISSUE}        # Jiho : NAMES
        ALLNAMES=${ALLNAMES}","${Sample_ID}"_"${TISSUE}
        CASE_BAM_PATH=${CASE_BAM_PATH}","${DATA_PATH}"/"${hg38}"/"${TISSUE}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
    done

    SAMTOOLS="/opt/Yonsei/samtools/1.7/"
    BCFTOOLS="/opt/Yonsei/bcftools/1.7/"
    BGZIP="~/miniconda3/envs/cnvpytor/bin/"
    SHAPEIT="~/miniconda3/envs/cnvpytor/bin/"
    PICARD="/opt/Yonsei/Picard/2.26.4/"
    MINCOV="50"
    MAXCOV="1000"
    READQUALITY="20"
    BASEQUALITY="20"
    REF_VERSION="hg38"
    RANDOM="1"
    dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.vcf.gz"
    CHR_NOTATION="true"
    PROCESSES=3
    MINREADS=50   #20 / 30 / 100 /20
    MAXREADS=1000  #1000 / 400 /3000 /600
    BIN="250kb"  #200kb
    PHASE="None"

    qsub -pe smp 5 -o $logPath"/hatchet_v1.0" -e $logPath"/hatchet_v1.0" -N "hat_v1.0_"$Sample_ID  -hold_jid "doc_"${Sample_ID}"_"${TISSUE} ${CURRENT_PATH}"/hatchet_v1.0_pipe.sh" \
        --CASE_BAM_PATH ${CASE_BAM_PATH} --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} --SAMPLES ${SAMPLES} --ALLNAMES ${ALLNAMES} \
        --HATCHET_PATH ${HATCHET_PATH} --REF ${REF_hg38} --REF_VERSION ${REF_VERSION} \
        --dbSNP ${dbSNP} --REGION ${INTERVAL} \
        --SAMTOOLS ${SAMTOOLS} --BCFTOOLS ${BCFTOOLS} --BGZIP ${BGZIP} --SHAPEIT ${SHAPEIT} --PICARD ${PICARD} \
        --MINCOV ${MINCOV} --MAXCOV ${MAXCOV} --MINREADS ${MINREADS} --MAXREADS ${MAXREADS}  --PROCESSES ${PROCESSES} --READQUALITY ${READQUALITY} --BASEQUALITY ${BASEQUALITY} \
        --CHR_NOTATION ${CHR_NOTATION} --BIN ${BIN} --PHASE ${PHASE} --RANDOM ${RANDOM} \
        --LOGPATH ${logPath} --SAMPLE_ID ${Sample_ID}

done


