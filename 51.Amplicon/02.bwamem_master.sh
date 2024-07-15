#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/51.Amplicon/02.multiplex"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in 11.bwa 12.postbwa 13.depthofcoverage; do
    # if [ $logPath"/"$sublog ] ; then
    #     rm -rf $logPath"/"$sublog
    # fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

################################################################# SINGLE ###################################################################    
DATA_PATH="/data/project/Meningioma/51.Amplicon/01.single"

for Sample_ID in 190426_NF2 220930_NF2 221026_NF2 230323_NF2 230920_NF2 230405_TRAF7; do  
    FQ=$(find "$DATA_PATH"/01.QC/00.raw/"${Sample_ID}" -name  "*1.fq.gz"  )
    FQ_LIST=(${FQ// / })                   # 이를 배열 (list)로 만듬

    if [ ! ${#FQ_LIST[@]} -eq 0 ]; then
        echo ${Sample_ID}
    fi

    for idx in ${!FQ_LIST[@]};do              # @ 배열의 모든 element    #! : indexing
        FASTQ_PATH_1=${FQ_LIST[idx]}         # idx번째의 파일명을 담아둔다
        FASTQ_PATH_2=${FASTQ_PATH_1%"1.fq.gz"}"2.fq.gz"

        S1=${FQ_LIST[idx]##*/}
        DATE_TISSUE=${S1%_1*.fq.gz}
        DATE=${S1%_[DVC]*}

        echo -e "\t"${DATE}"\t"${DATE_TISSUE}

        FASTP_PATH_1=$DATA_PATH"/01.QC/01.fastp/"${Sample_ID}"/"${DATE_TISSUE}".R1.fq.gz"
        FASTP_PATH_2=$DATA_PATH"/01.QC/01.fastp/"${Sample_ID}"/"${DATE_TISSUE}".R2.fq.gz"
        FASTQC_PATH=$DATA_PATH"/01.QC/02.fastqc/"${Sample_ID}

       
        REF_hg38="/home/goldpm1/reference/genome.fa"
        #REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
        REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"

        if [ "${Sample_ID}" == "230405_TRAF7" ]; then
            REF="/home/goldpm1/reference/chr16/chr16.fa"
            #REF="/home/goldpm1/reference/TRAF7/TRAF7.genome.fa"
        else
            REF="${REF_hg38}"
        fi

        hg="hg38"

        PRE_BAM_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${Sample_ID}"/01.Pre_bam/"${DATE_TISSUE}".sorted.bam"
        MarkDuplicate_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${Sample_ID}"/02.MarkDuplicate/"${DATE_TISSUE}".mkdp.sorted.bam"
        AddOrReplaceReadGroups_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${Sample_ID}"/03.AddOrReplaceReadGroups/"${DATE_TISSUE}".arg.mkdp.sorted.bam"
        BQSR_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${Sample_ID}"/04.BQSR/"${DATE_TISSUE}".recal.arg.mkdp.sorted.bam"
        BQSR_RECAL_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${Sample_ID}"/04.BQSR/"${DATE_TISSUE}".time.recal.table"
        LOCAL_REALIGNMENT_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${Sample_ID}"/05.LocalRealignment/"${DATE_TISSUE}".bam"
        FINAL_BAM_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${Sample_ID}"/06.Final_bam/"${DATE_TISSUE}".bam"
        DOC_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${Sample_ID}"/07.DepthOfCoverage/"${DATE_TISSUE}".depth.result"
        

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
        
        
        qsub -pe smp 5 -e ${logPath}"/11.bwa" -o ${logPath}"/11.bwa" -N "bwa_"${DATE_TISSUE} -hold_jid 'FP_'${DATE_TISSUE}",FQC_"${DATE_TISSUE} "02.bwamem_pipe_01.bwa.sh" \
        --FASTQ_PATH_1 ${FASTP_PATH_1} --FASTQ_PATH_2 ${FASTP_PATH_2}  --PRE_BAM_PATH ${PRE_BAM_PATH} --REF ${REF} --TMP_PATH ${TMP_PATH}
        
        qsub -pe smp 8 -e ${logPath}"/12.postbwa" -o ${logPath}"/12.postbwa" -N "postbwa_"${DATE_TISSUE} -hold_jid 'bwa_'${DATE_TISSUE} "02.bwamem_pipe_02.postbwa.sh" \
        --DATE_TISSUE ${DATE_TISSUE} --PRE_BAM_PATH ${PRE_BAM_PATH} --MarkDuplicate_PATH ${MarkDuplicate_PATH} --AddOrReplaceReadGroups_PATH ${AddOrReplaceReadGroups_PATH} \
        --BQSR_PATH ${BQSR_PATH} --BQSR_RECAL_PATH ${BQSR_RECAL_PATH}  --LOCAL_REALIGNMENT_PATH ${LOCAL_REALIGNMENT_PATH} --FINAL_BAM_PATH ${FINAL_BAM_PATH} \
        --REF ${REF} --dbsnp ${dbsnp} --TMP_PATH ${TMP_PATH}
        
        # qsub -pe smp 6 -e ${logPath}"/13.depthofcoverage" -o ${logPath}"/13.depthofcoverage" -N "doc_"${DATE_TISSUE} -hold_jid 'postbwa_'${DATE_TISSUE} "02.bwamem_pipe_03.depthofcoverage.sh" \
        # --REF ${REF_hg38}  --DOC_PATH ${DOC_PATH} --FINAL_BAM_PATH ${FINAL_BAM_PATH}  --INTERVAL ${INTERVAL}
    
    done    
done





################################################################# MULTIPLEX ###################################################################    

# DATA_PATH="/data/project/Meningioma/51.Amplicon/02.multiplex"

# for Sample_ID in 221026_Dura; do
#     FASTQ_PATH_1=$(find "$DATA_PATH"/01.QC/00.raw/ -name  ${Sample_ID}*1.fq.gz  )
#     FASTQ_PATH_2=${FASTQ_PATH_1%"1.fq.gz"}"2.fq.gz"

#     S1=${FASTQ_PATH_1##*/}
#     DATE_TISSUE=${S1%_1*.fq.gz}
#     DATE=${S1%_[DVC]*}

    
#     echo -e "Multiplex\t"${DATE_TISSUE}

#     FASTP_PATH_1=$DATA_PATH"/01.QC/01.fastp/"${DATE_TISSUE}".R1.fq.gz"
#     FASTP_PATH_2=$DATA_PATH"/01.QC/01.fastp/"${DATE_TISSUE}".R2.fq.gz"

#     REF_hg38="/home/goldpm1/reference/genome.fa"
#     #REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
#     REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"

#     hg="hg38"

#     PRE_BAM_PATH=${DATA_PATH}"/02.Align/"${hg}"//01.Pre_bam/"${DATE_TISSUE}".sorted.bam"
#     MarkDuplicate_PATH=${DATA_PATH}"/02.Align/"${hg}"//02.MarkDuplicate/"${DATE_TISSUE}".mkdp.sorted.bam"
#     AddOrReplaceReadGroups_PATH=${DATA_PATH}"/02.Align/"${hg}"//03.AddOrReplaceReadGroups/"${DATE_TISSUE}".arg.mkdp.sorted.bam"
#     BQSR_PATH=${DATA_PATH}"/02.Align/"${hg}"//04.BQSR/"${DATE_TISSUE}".recal.arg.mkdp.sorted.bam"
#     BQSR_RECAL_PATH=${DATA_PATH}"/02.Align/"${hg}"//04.BQSR/"${DATE_TISSUE}".time.recal.table"
#     LOCAL_REALIGNMENT_PATH=${DATA_PATH}"/02.Align/"${hg}"//05.LocalRealignment/"${DATE_TISSUE}".bam"
#     FINAL_BAM_PATH=${DATA_PATH}"/02.Align/"${hg}"//06.Final_bam/"${DATE_TISSUE}".bam"
#     DOC_PATH=${DATA_PATH}"/02.Align/"${hg}"//07.DepthOfCoverage/"${DATE_TISSUE}".depth.result"
    

#     dbsnp="/data/public/dbSNP/b155/GRCh38/GCF_000001405.39.re.vcf.gz"
#     TMP_PATH=${DATA_PATH}"/temp"
#     INTERVAL="/home/goldpm1/resources/Exon.reference.GRCh38.bed"
#     INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
#     INTERVAL="/home/goldpm1/resources/TMB359.theragen.hg38.bed"
    
#     for folder in ${PRE_BAM_PATH%/*} ${MarkDuplicate_PATH%/*} ${AddOrReplaceReadGroups_PATH%/*} ${BQSR_PATH%/*} ${BQSR_RECAL_PATH%/*} ${LOCAL_REALIGNMENT_PATH%/*} ${FINAL_BAM_PATH%/*} ${DOC_PATH%/*} ${TMP_PATH}; do
#         if [ ! -d $folder ] ; then
#             mkdir -p $folder
#         fi
#     done
    

#     qsub -pe smp 5 -e ${logPath}"/11.bwa" -o ${logPath}"/11.bwa" -N "bwa_"${DATE_TISSUE} -hold_jid 'FP_'${DATE_TISSUE}",FQC_"${DATE_TISSUE} "02.bwamem_pipe_01.bwa.sh" \
#     --FASTQ_PATH_1 ${FASTP_PATH_1} --FASTQ_PATH_2 ${FASTP_PATH_2}  --PRE_BAM_PATH ${PRE_BAM_PATH} --REF ${REF_hg38} --TMP_PATH ${TMP_PATH}
    
#     qsub -pe smp 8 -e ${logPath}"/12.postbwa" -o ${logPath}"/12.postbwa" -N "postbwa_"${DATE_TISSUE} -hold_jid 'bwa_'${DATE_TISSUE} "02.bwamem_pipe_02.postbwa.sh" \
#     --DATE_TISSUE ${DATE_TISSUE} --PRE_BAM_PATH ${PRE_BAM_PATH} --MarkDuplicate_PATH ${MarkDuplicate_PATH} --AddOrReplaceReadGroups_PATH ${AddOrReplaceReadGroups_PATH} \
#     --BQSR_PATH ${BQSR_PATH} --BQSR_RECAL_PATH ${BQSR_RECAL_PATH}  --LOCAL_REALIGNMENT_PATH ${LOCAL_REALIGNMENT_PATH} --FINAL_BAM_PATH ${FINAL_BAM_PATH} \
#     --REF ${REF_hg38} --dbsnp ${dbsnp} --TMP_PATH ${TMP_PATH}
    
#     # qsub -pe smp 6 -e ${logPath}"/13.depthofcoverage" -o ${logPath}"/13.depthofcoverage" -N "doc_"${DATE_TISSUE} -hold_jid 'postbwa_'${DATE_TISSUE} "02.bwamem_pipe_03.depthofcoverage.sh" \
#     # --REF ${REF_hg38}  --DOC_PATH ${DOC_PATH} --FINAL_BAM_PATH ${FINAL_BAM_PATH}  --INTERVAL ${INTERVAL}
    

# done

