#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/00.Targetted"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in 11.bwa 12.postbwa 13.depthofcoverage; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done



# sample_name_list=$(cat ${CURRENT_PATH}"/sample_name.txt")
# sample_name_LIST=(${sample_name_list// / })     # array로 만듬
# for idx in ${!sample_name_LIST[@]}; do
#     Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102
    
# Align
for Sample_ID in 230526_2 230526_3; do
    for TISSUE in Dura; do
        FASTQ_PATH_1=${DATA_PATH}"/01.QC/01.fastp/"${TISSUE}"/"${Sample_ID}"_"${TISSUE}".R1.fq.gz"
        FASTQ_PATH_2=${DATA_PATH}"/01.QC/01.fastp/"${TISSUE}"/"${Sample_ID}"_"${TISSUE}".R2.fq.gz"
        
        REF_hg38="/home/goldpm1/reference/genome.fa"
        #REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
        REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"

        hg="hg38"

        PRE_BAM_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${TISSUE}"/01.Pre_bam/"${Sample_ID}".sorted.bam"
        MarkDuplicate_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${TISSUE}"/02.MarkDuplicate/"${Sample_ID}".mkdp.sorted.bam"
        AddOrReplaceReadGroups_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${TISSUE}"/03.AddOrReplaceReadGroups/"${Sample_ID}".arg.mkdp.sorted.bam"
        BQSR_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${TISSUE}"/04.BQSR/"${Sample_ID}".recal.arg.mkdp.sorted.bam"
        BQSR_RECAL_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${TISSUE}"/04.BQSR/"${Sample_ID}".time.recal.table"
        LOCAL_REALIGNMENT_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${TISSUE}"/05.LocalRealignment/"${Sample_ID}"_"${TISSUE}".bam"
        FINAL_BAM_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${TISSUE}"/06.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
        DOC_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${TISSUE}"/07.DepthOfCoverage/"${Sample_ID}"_"${TISSUE}".depth.result"
        

        dbsnp="/data/public/dbSNP/b155/GRCh38/GCF_000001405.39.re.vcf.gz"
        TMP_PATH=${DATA_PATH}"/temp"
        INTERVAL="/home/goldpm1/resources/Exon.reference.GRCh38.bed"
        INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
        INTERVAL="/home/goldpm1/resources/TMB359.theragen.hg38.bed"
        
        for folder in ${PRE_BAM_PATH%/*} ${MarkDuplicate_PATH%/*} ${AddOrReplaceReadGroups_PATH%/*} ${BQSR_PATH%/*} ${BQSR_RECAL_PATH%/*} ${LOCAL_REALIGNMENT_PATH%/*} ${FINAL_BAM_PATH%/*} ${DOC_PATH%/*} ${TMP_PATH}; do
            if [ ! -d $folder ] ; then
                mkdir -p $folder
            fi
        done
        
        
        qsub -pe smp 5 -e ${logPath}"/11.bwa" -o ${logPath}"/11.bwa" -N "bwa_"${Sample_ID}"_"${TISSUE} -hold_jid 'FP_'${Sample_ID}"_"${TISSUE}",FQC_"${Sample_ID}"_"${TISSUE} "02.bwamem_pipe_01.bwa.sh" \
        --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2}  --PRE_BAM_PATH ${PRE_BAM_PATH} --REF ${REF_hg38} --TMP_PATH ${TMP_PATH}
        
        qsub -pe smp 8 -e ${logPath}"/12.postbwa" -o ${logPath}"/12.postbwa" -N "postbwa_"${Sample_ID}"_"${TISSUE} -hold_jid 'bwa_'${Sample_ID}"_"${TISSUE} "02.bwamem_pipe_02.postbwa.sh" \
        --Sample_ID ${Sample_ID}  --TISSUE ${TISSUE} --PRE_BAM_PATH ${PRE_BAM_PATH} --MarkDuplicate_PATH ${MarkDuplicate_PATH} --AddOrReplaceReadGroups_PATH ${AddOrReplaceReadGroups_PATH} \
        --BQSR_PATH ${BQSR_PATH} --BQSR_RECAL_PATH ${BQSR_RECAL_PATH}  --LOCAL_REALIGNMENT_PATH ${LOCAL_REALIGNMENT_PATH} --FINAL_BAM_PATH ${FINAL_BAM_PATH} \
        --REF ${REF_hg38} --dbsnp ${dbsnp} --TMP_PATH ${TMP_PATH}
        
        qsub -pe smp 6 -e ${logPath}"/13.depthofcoverage" -o ${logPath}"/13.depthofcoverage" -N "doc_"${Sample_ID}"_"${TISSUE} -hold_jid 'postbwa_'${Sample_ID}"_"${TISSUE} "02.bwamem_pipe_03.depthofcoverage.sh" \
        --REF ${REF_hg38}  --DOC_PATH ${DOC_PATH} --FINAL_BAM_PATH ${FINAL_BAM_PATH}  --INTERVAL ${INTERVAL}
        
    done
done
