#!/bin/bash
#$ -cwd
#$ -S /bin/bash


############### Somatic matched pair ##################


REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
REF_hg19="/home/goldpm1/reference/Broadhg19/Homo_sapiens_assembly19.fasta"             # 왜그런지는 모르겠지만 이 reference를 써야 돌아간다
#REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
REF_hg38="/home/goldpm1/reference/genome.fa"

INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed.gz"
dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz"


CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

for sublog in "01.manta"; do
    if [ -d $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

DATA_PATH="/data/project/Meningioma/02.Align"
MANTA_PATH="/data/project/Meningioma/21.SV/01.manta"

for subpath in "01.raw" "02.PASS" "03.pandas"; do
    if [ ! -d $MANTA_PATH"/"$subpath ] ; then
        mkdir -p $MANTA_PATH"/"$subpath
    fi
done

sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬


for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102

    CONTROL_BAM_PATH=${DATA_PATH}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"

    for TISSUE in Tumor Dura ; do   #Tumor Dura
        CASE_BAM_PATH=${DATA_PATH}"/"${TISSUE}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"

        OUTPUT_DIR=${MANTA_PATH}"/01.raw/"${TISSUE}"/"${Sample_ID}
        OUTPUT_PASS_PATH=${MANTA_PATH}"/02.PASS/"${TISSUE}"/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".Manta.PASS.vcf"
        if [ -d ${OUTPUT_DIR} ] ; then
            rm -rf ${OUTPUT_DIR}
        fi
        if [ ! -d ${OUTPUT_DIR} ] ; then
            mkdir -p ${OUTPUT_DIR}
        fi
        if [ -d ${OUTPUT_PASS_PATH%/*} ] ; then
            rm -rf ${OUTPUT_PASS_PATH%/*}
        fi
        if [ ! -d ${OUTPUT_PASS_PATH%/*} ] ; then
            mkdir -p ${OUTPUT_PASS_PATH%/*}
        fi

        #01. Manta configuration & Execution
        # qsub -pe smp 3 -o $logPath"/01.manta" -e $logPath"/01.manta" -N "man01_"${Sample_ID}"_"${TISSUE}  "manta_pipe_01.sh" \
        #     --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} --CASE_BAM_PATH ${CASE_BAM_PATH}   \
        #     --REF ${REF_hg38} --CALLREGIONS ${INTERVAL} --OUTPUT_DIR ${OUTPUT_DIR} --OUTPUT_PASS_PATH ${OUTPUT_PASS_PATH}



        # 02. Pandas Dataframe으로 만들기

        INPUT_PATH=${MANTA_PATH}"/02.PASS/"${TISSUE}"/"${Sample_ID}"_"${TISSUE}".Manta.PASS.vcf"
        OUTPUT_PATH=${MANTA_PATH}"/02.PASS/"${TISSUE}"/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".Manta.PASS.chr.vcf"
        OUTPUT_INV_PATH=${MANTA_PATH}"/02.PASS/"${TISSUE}"/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".Manta.PASS.INV.vcf"
        PANDAS_DIR=${MANTA_PATH}"/03.pandas/"${TISSUE}"/"${Sample_ID}
        if [ -d ${PANDAS_DIR} ] ; then
            rm -rf ${PANDAS_DIR}
        fi
        if [ ! -d ${PANDAS_DIR} ] ; then
            mkdir -p ${PANDAS_DIR}
        fi
        # qsub -pe smp 1 -o $logPath"/01.manta" -e $logPath"/01.manta" -N "man02_"${Sample_ID}"_"${TISSUE}  "manta_pipe_02.sh" \
        #     --INPUT_PATH ${INPUT_PATH} --OUTPUT_PATH ${OUTPUT_PATH} --OUTPUT_INV_PATH ${OUTPUT_INV_PATH} --PANDAS_DIR ${PANDAS_DIR} \
        #     --TISSUE ${TISSUE} --Sample_ID ${Sample_ID} --REF ${REF_hg38} 

    done
done


