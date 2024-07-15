#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/02.Align"

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

REF_hg38="/home/goldpm1/reference/genome.fa"
#REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"

hg="hg38"
REF=${REF_hg38}

# sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
# sample_name_LIST=(${sample_name_list// / })     # array로 만듬
# for idx in ${!sample_name_LIST[@]}; do
#     Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102


#for Sample_ID in 220930 221026 221102 221202 230127 230323_2 230405_2 230419 230526 230822 230920; do
for Sample_ID in 230303; do 
    
    for TISSUE in Dura; do
    #for TISSUE in Blood Tumor  Tumor_FFT Tumor_PCT Tumor_PFT Tumor_PP Tumor_PT Dura Ventricle Cortex; do  
        FASTQ_PATH_1=${DATA_PATH%/*}"/01.QC/01.fastp/"${TISSUE%_*}"/"${Sample_ID}"_"${TISSUE}".R1.fq.gz"
        FASTQ_PATH_2=${DATA_PATH%/*}"/01.QC/01.fastp/"${TISSUE%_*}"/"${Sample_ID}"_"${TISSUE}".R2.fq.gz"
        
        if [ -f "$FASTQ_PATH_1" ]; then     # File이 있어야만 진행
            PRE_BAM_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE%_*}"/01.Pre_bam/"${Sample_ID}"_"${TISSUE}".sorted.bam"
            MarkDuplicate_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE%_*}"/02.MarkDuplicate/"${Sample_ID}"_"${TISSUE}".mkdp.sorted.bam"
            AddOrReplaceReadGroups_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE%_*}"/03.AddOrReplaceReadGroups/"${Sample_ID}"_"${TISSUE}".arg.mkdp.sorted.bam"
            BQSR_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE%_*}"/04.BQSR/"${Sample_ID}"_"${TISSUE}".recal.arg.mkdp.sorted.bam"
            BQSR_RECAL_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE%_*}"/04.BQSR/"${Sample_ID}"_"${TISSUE}".time.recal.table"
            FINAL_BAM_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE%_*}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
            DOC_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE%_*}"/06.DepthOfCoverage/"${Sample_ID}"_"${TISSUE}".depth.result"
            

            dbsnp="/data/public/dbSNP/b155/GRCh38/GCF_000001405.39.re.vcf.gz"
            TMP_PATH=${DATA_PATH}"/temp"
            INTERVAL="/home/goldpm1/resources/Exon.reference.GRCh38.bed"
            INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
            
            for folder in ${PRE_BAM_PATH%/*} ${MarkDuplicate_PATH%/*} ${AddOrReplaceReadGroups_PATH%/*} ${BQSR_PATH%/*} ${BQSR_RECAL_PATH%/*} ${FINAL_BAM_PATH%/*} ${DOC_PATH%/*} ${TMP_PATH%/*}; do
                if [ ! -d $folder ] ; then
                    mkdir -p $folder
                fi
            done
            
            
            qsub -pe smp 5 -e ${logPath}"/01.bwa" -o ${logPath}"/01.bwa" -N "bwa_"${Sample_ID}"_"${TISSUE} -hold_jid 'FP_'${Sample_ID}"_"${TISSUE}",FQC_"${Sample_ID}"_"${TISSUE} bwamem_pipe_01.bwa.sh \
            --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2}  --PRE_BAM_PATH ${PRE_BAM_PATH} --REF ${REF} --TMP_PATH ${TMP_PATH}
        
            qsub -pe smp 8 -e ${logPath}"/02.postbwa" -o ${logPath}"/02.postbwa" -N "postbwa_"${Sample_ID}"_"${TISSUE} -hold_jid 'bwa_'${Sample_ID}"_"${TISSUE} bwamem_pipe_02.postbwa.sh \
            --Sample_ID ${Sample_ID}  --TISSUE ${TISSUE} --PRE_BAM_PATH ${PRE_BAM_PATH} --MarkDuplicate_PATH ${MarkDuplicate_PATH} --AddOrReplaceReadGroups_PATH ${AddOrReplaceReadGroups_PATH} \
            --BQSR_PATH ${BQSR_PATH} --BQSR_RECAL_PATH ${BQSR_RECAL_PATH}  --FINAL_BAM_PATH ${FINAL_BAM_PATH} \
            --REF ${REF} --dbsnp ${dbsnp} --TMP_PATH ${TMP_PATH}
            
            # qsub -pe smp 5 -e ${logPath}"/03.depthofcoverage" -o ${logPath}"/03.depthofcoverage" -N "doc_"${Sample_ID}"_"${TISSUE} -hold_jid 'postbwa_'${Sample_ID}"_"${TISSUE} depthofcoverage.sh \
            # --REF ${REF}  --DOC_PATH ${DOC_PATH} --FINAL_BAM_PATH ${FINAL_BAM_PATH}  --INTERVAL ${INTERVAL}
        fi    
    done
done
