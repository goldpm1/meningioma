#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/61.Lowinput/01.XT_HS"
#DATA_PATH="/data/project/Meningioma/61.Lowinput/02.PTA"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "02-1.bwa" "02-2.postbwa" "02-3.depthofcoverage"; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done



for Sample_ID in 190426; do
    for Clone_No in 3 4; do
        FASTQ_PATH_1=${DATA_PATH}"/01.QC/01.fastp/Tumor/"${Sample_ID}"_"${Clone_No}".R1.fq.gz"
        FASTQ_PATH_2=${DATA_PATH}"/01.QC/01.fastp/Tumor/"${Sample_ID}"_"${Clone_No}".R2.fq.gz"

        if [ ! -f ${FASTQ_PATH_1} ] ; then
            continue  # Skip to the next iteration
        fi
        echo ${Sample_ID}"_"${Clone_No}
        
        REF_hg38="/home/goldpm1/reference/genome.fa"
        #REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
        REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"

        hg="hg38"

        PRE_BAM_PATH=${DATA_PATH}"/02.Align/"${hg}"/Tumor/01.Pre_bam/"${Sample_ID}"_"${Clone_No}".sorted.bam"
        MarkDuplicate_PATH=${DATA_PATH}"/02.Align/"${hg}"/Tumor/02.MarkDuplicate/"${Sample_ID}"_"${Clone_No}".mkdp.sorted.bam"
        AddOrReplaceReadGroups_PATH=${DATA_PATH}"/02.Align/"${hg}"/Tumor/03.AddOrReplaceReadGroups/"${Sample_ID}"_"${Clone_No}".arg.mkdp.sorted.bam"
        BQSR_PATH=${DATA_PATH}"/02.Align/"${hg}"/Tumor/04.BQSR/"${Sample_ID}"_"${Clone_No}".recal.arg.mkdp.sorted.bam"
        BQSR_RECAL_PATH=${DATA_PATH}"/02.Align/"${hg}"/Tumor/04.BQSR/"${Sample_ID}"_"${Clone_No}".time.recal.table"
        LOCAL_REALIGNMENT_PATH=${DATA_PATH}"/02.Align/"${hg}"/Tumor/05.LocalRealignment/"${Sample_ID}"_"${Clone_No}".bam"
        FINAL_BAM_PATH=${DATA_PATH}"/02.Align/"${hg}"/Tumor/06.Final_bam/"${Sample_ID}"_"${Clone_No}".bam"
        DOC_PATH=${DATA_PATH}"/02.Align/"${hg}"/Tumor/07.DepthOfCoverage/"${Sample_ID}"_"${Clone_No}".depth.result"
        

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
        
        
        qsub -pe smp 5 -e ${logPath}"/02-1.bwa" -o ${logPath}"/02-1.bwa" -N "bwa_"${Sample_ID}"_"${Clone_No} -hold_jid 'FP_'${Sample_ID}"_"${Clone_No}",FQC_"${Sample_ID}"_"${Clone_No} "02.bwamem_pipe_01.bwa.sh" \
        --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2}  --PRE_BAM_PATH ${PRE_BAM_PATH} --REF ${REF_hg38} --TMP_PATH ${TMP_PATH}
        
        qsub -pe smp 8 -e ${logPath}"/02-2.postbwa" -o ${logPath}"/02-2.postbwa" -N "postbwa_"${Sample_ID}"_"${Clone_No} -hold_jid 'bwa_'${Sample_ID}"_"${Clone_No} "02.bwamem_pipe_02.postbwa.sh" \
        --Sample_ID ${Sample_ID}  --TISSUE ${Clone_No} --PRE_BAM_PATH ${PRE_BAM_PATH} --MarkDuplicate_PATH ${MarkDuplicate_PATH} --AddOrReplaceReadGroups_PATH ${AddOrReplaceReadGroups_PATH} \
        --BQSR_PATH ${BQSR_PATH} --BQSR_RECAL_PATH ${BQSR_RECAL_PATH}  --LOCAL_REALIGNMENT_PATH ${LOCAL_REALIGNMENT_PATH} --FINAL_BAM_PATH ${FINAL_BAM_PATH} \
        --REF ${REF_hg38} --dbsnp ${dbsnp} --TMP_PATH ${TMP_PATH}
        
        qsub -pe smp 6 -e ${logPath}"/02-3.depthofcoverage" -o ${logPath}"/02-3.depthofcoverage" -N "doc_"${Sample_ID}"_"${Clone_No} -hold_jid 'postbwa_'${Sample_ID}"_"${Clone_No} "02.bwamem_pipe_03.depthofcoverage.sh" \
        --REF ${REF_hg38}  --DOC_PATH ${DOC_PATH} --FINAL_BAM_PATH ${FINAL_BAM_PATH}  --INTERVAL ${INTERVAL}
        
    done
done
