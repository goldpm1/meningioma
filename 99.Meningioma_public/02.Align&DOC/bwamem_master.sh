#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

PROJECT_DIR="/data/project/Meningioma/99.Meningioma_public/SRP261564"     # SRP050339 (single)  SRP227246 (Al-mefty) SRP261564 (Saudi, single)
FASTQ_DIR=${PROJECT_DIR}"/02.fastq"
BAM_DIR=${PROJECT_DIR}"/03.bam"

hg="hg38"
if [ "${hg}" == "hg38" ]; then
    REF="/home/goldpm1/reference/genome.fa"           #REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
    dbsnp="/data/public/dbSNP/b155/GRCh38/GCF_000001405.39.re.vcf.gz"
    INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
elif [ "${hg}" == "hg19" ]; then
    REF="/home/goldpm1/reference/hg19/hg19.fa"
    dbsnp="/data/public/dbSNP/b155/GRCh37/GCF_000001405.25.re.vcf.gz"
    INTERVAL="/home/goldpm1/resources/Exon.reference.GRCh38.bed"
fi



if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in 01.bwa 02.postbwa 03.depthofcoverage; do
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
    
    # FASTQ_PATH_1=${FASTQ_DIR}"/"${Sample_ID}".R1.fq.gz"
    # FASTQ_PATH_2=${FASTQ_DIR}"/"${Sample_ID}".R2.fq.gz"

    FASTQ_PATH_1=${FASTQ_DIR}"/"${Sample_ID}".fastq"
    
    if [ -f "$FASTQ_PATH_1" ]; then     # File이 있어야만 진행
        PRE_BAM_PATH=${BAM_DIR}"/"${hg}"/01.Pre_bam/"${Sample_ID}".sorted.bam"
        MarkDuplicate_PATH=${BAM_DIR}"/"${hg}"/02.MarkDuplicate/"${Sample_ID}".mkdp.sorted.bam"
        AddOrReplaceReadGroups_PATH=${BAM_DIR}"/"${hg}"/03.AddOrReplaceReadGroups/"${Sample_ID}".arg.mkdp.sorted.bam"
        BQSR_PATH=${BAM_DIR}"/"${hg}"/04.BQSR/"${Sample_ID}".recal.arg.mkdp.sorted.bam"
        BQSR_RECAL_PATH=${BAM_DIR}"/"${hg}"/04.BQSR/"${Sample_ID}".time.recal.table"
        FINAL_BAM_PATH=${BAM_DIR}"/"${hg}"/05.Final_bam/"${Sample_ID}".bam"
        DOC_PATH=${BAM_DIR}"/"${hg}"/06.DepthOfCoverage/"${Sample_ID}".depth.result"
        TMP_PATH=${BAM_DIR}"/temp"

        
        for folder in ${PRE_BAM_PATH%/*} ${MarkDuplicate_PATH%/*} ${AddOrReplaceReadGroups_PATH%/*} ${BQSR_PATH%/*} ${BQSR_RECAL_PATH%/*} ${FINAL_BAM_PATH%/*} ${DOC_PATH%/*} ${TMP_PATH%/*}; do
            if [ ! -d $folder ] ; then
                mkdir -p $folder
            fi
        done
        
        
        # qsub -pe smp 5 -e ${logPath}"/01.bwa" -o ${logPath}"/01.bwa" -N "bwa_"${Sample_ID} -hold_jid 'FP_'${Sample_ID}",FQC_"${Sample_ID} bwamem_pipe_01.bwa.sh \
        # --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2}  --PRE_BAM_PATH ${PRE_BAM_PATH} --REF ${REF} --TMP_PATH ${TMP_PATH}

        # single end
        # qsub -pe smp 5 -e ${logPath}"/01.bwa" -o ${logPath}"/01.bwa" -N "bwa_"${Sample_ID} -hold_jid 'FP_'${Sample_ID}",FQC_"${Sample_ID} bwamem_pipe_01.bwa_singleend.sh \
        # --FASTQ_PATH_1 ${FASTQ_PATH_1}   --PRE_BAM_PATH ${PRE_BAM_PATH} --REF ${REF} --TMP_PATH ${TMP_PATH}
        
        qsub -pe smp 8 -e ${logPath}"/02.postbwa" -o ${logPath}"/02.postbwa" -N "postbwa_"${Sample_ID} -hold_jid 'bwa_'${Sample_ID} bwamem_pipe_02.postbwa.sh \
        --Sample_ID ${Sample_ID}   --PRE_BAM_PATH ${PRE_BAM_PATH} --MarkDuplicate_PATH ${MarkDuplicate_PATH} --AddOrReplaceReadGroups_PATH ${AddOrReplaceReadGroups_PATH} \
        --BQSR_PATH ${BQSR_PATH} --BQSR_RECAL_PATH ${BQSR_RECAL_PATH}  --FINAL_BAM_PATH ${FINAL_BAM_PATH} \
        --REF ${REF} --dbsnp ${dbsnp} --TMP_PATH ${TMP_PATH}
        
        # qsub -pe smp 5 -e ${logPath}"/03.depthofcoverage" -o ${logPath}"/03.depthofcoverage" -N "doc_"${Sample_ID} -hold_jid 'postbwa_'${Sample_ID} depthofcoverage.sh \
        # --REF ${REF}  --DOC_PATH ${DOC_PATH} --FINAL_BAM_PATH ${FINAL_BAM_PATH}  --INTERVAL ${INTERVAL}
    fi    
    
done
