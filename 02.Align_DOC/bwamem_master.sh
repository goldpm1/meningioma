#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/02.Align"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in 01.bwa 02.postbwa 03.depthofcoverage 04.mosdepth; do
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



# Tumor only
#for Sample_ID in 250310 250319 250407 250408 250428 250509 250513 250515 250520 250526_Left 250526_Right 250527 250530  250602 250605 250609  \
# 221021 231011 231208 231220 240112 250227 250617 250627 250704 250716 250723_Ant 250723_Post 250725; do 

# Paired
# for Sample_ID in 190426 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT 241127 241127_FT 241127_FTRt 241127_Para \
#     190524 220930_2 221026 221102 221202 230127 230303_1 230323_2 230323_11 230405_2 230419 230526 230812 230822 230920 231006 231025 231101 231206 240110 240320 240325 240612 240911 241023 \
#     230802_2 240403_3 240412_2; do 

# for Sample_ID in  241016 241211 250425 250212 250502 250509 \
#     250319 250530 250408 250513 250520 \
#     220930_3 221021 231011 231208 231220 240112 250227 250617 250627 250704 250716 250723_Ant 250723_Post 250725; do 
# for Sample_ID in  220930 230405 230303 230323 230802 240403 240412; do 

for Sample_ID in SMS_1st SMS_2nd SMS_3rd SMS_4th SMS_mets SMS; do 
    echo $Sample_ID
    for TISSUE in Tumor Dura Falx DT Falx_1 Falx_2 Ventricle Cortex Blood; do     # Ventricle Cortex
    #for TISSUE in Tumor Dura Dura_Ant Blood; do     # Falx_1 Falx_2 Falx_FT Ventricle Cortex 
        if [[ $TISSUE == *"Falx"* || $TISSUE == *"Falx_1"* || $TISSUE == *"Falx_2"* || $TISSUE == *"Falx_FT"* || $TISSUE == *"Dura"* || $TISSUE == *"DT"* ]]; then
            TISSUE_BROAD="Dura"
        elif  [[ $TISSUE == "Tumor" ]]; then
            TISSUE_BROAD="Tumor"
        else
            TISSUE_BROAD=${TISSUE%%_*}
        fi

        FASTQ_PATH_1=${DATA_PATH%/*}"/01.QC/01.fastp/"${TISSUE_BROAD}"/"${Sample_ID}"_"${TISSUE}".R1.fq.gz"
        FASTQ_PATH_2=${DATA_PATH%/*}"/01.QC/01.fastp/"${TISSUE_BROAD}"/"${Sample_ID}"_"${TISSUE}".R2.fq.gz"
        
        if [ -f "$FASTQ_PATH_1" ]; then     # File이 있어야만 진행
            PRE_BAM_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE_BROAD}"/01.Pre_bam/"${Sample_ID}"_"${TISSUE}".sorted.bam"
            MarkDuplicate_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE_BROAD}"/02.MarkDuplicate/"${Sample_ID}"_"${TISSUE}".mkdp.sorted.bam"
            AddOrReplaceReadGroups_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE_BROAD}"/03.AddOrReplaceReadGroups/"${Sample_ID}"_"${TISSUE}".arg.mkdp.sorted.bam"
            BQSR_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE_BROAD}"/04.BQSR/"${Sample_ID}"_"${TISSUE}".recal.arg.mkdp.sorted.bam"
            BQSR_RECAL_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE_BROAD}"/04.BQSR/"${Sample_ID}"_"${TISSUE}".time.recal.table"
            FINAL_BAM_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE_BROAD}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
            DOC_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE_BROAD}"/06.DepthOfCoverage/"${Sample_ID}"_"${TISSUE}".depth.result"
            MOSDEPTH_DIR=${DATA_PATH}"/"${hg}"/"${TISSUE_BROAD}"/07.mosdepth"
            

            dbsnp="/data/public/dbSNP/b155/GRCh38/GCF_000001405.39.re.vcf.gz"
            TMP_PATH=${DATA_PATH}"/temp"
            INTERVAL="/home/goldpm1/resources/Exon.reference.GRCh38.bed"
            INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
            
            for folder in ${PRE_BAM_PATH%/*} ${MarkDuplicate_PATH%/*} ${AddOrReplaceReadGroups_PATH%/*} ${BQSR_PATH%/*} ${BQSR_RECAL_PATH%/*} ${FINAL_BAM_PATH%/*} ${DOC_PATH%/*} ${TMP_PATH%/*} ${MOSDEPTH_DIR}; do
                if [ ! -d $folder ] ; then
                    mkdir -p $folder
                fi
            done
            
            
            BWA_JOB=$(sbatch --job-name='bwa_'${Sample_ID}"_"${TISSUE}"_"${hg} --output=${logPath}"/01.bwa/"${Sample_ID}"_"${TISSUE}".out" --error=${logPath}"/01.bwa/"${Sample_ID}"_"${TISSUE}".err" \
                --cpus-per-task=5 --mem-per-cpu=14GB --time=4:00:00 --partition=cpu --qos nstumor \
                --wrap="bash bwamem_pipe_01.bwa.sh --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2}  --PRE_BAM_PATH ${PRE_BAM_PATH} --REF ${REF} --TMP_PATH ${TMP_PATH}")
            BWA_JOBID=$(echo ${BWA_JOB} | awk '{print $NF}')

            POSTBWA_JOB=$(sbatch --dependency=afterok:${BWA_JOBID} --job-name='postbwa_'${Sample_ID}"_"${TISSUE}"_"${hg} --output=${logPath}"/02.postbwa/"${Sample_ID}"_"${TISSUE}".out" --error=${logPath}"/02.postbwa/"${Sample_ID}"_"${TISSUE}".err" \
                --cpus-per-task=8 --mem-per-cpu=14GB --time=10:00:00 --partition=cpu --qos nstumor \
                --wrap="bash bwamem_pipe_02.postbwa.sh --Sample_ID ${Sample_ID}  --TISSUE ${TISSUE} --PRE_BAM_PATH ${PRE_BAM_PATH} --MarkDuplicate_PATH ${MarkDuplicate_PATH} --AddOrReplaceReadGroups_PATH ${AddOrReplaceReadGroups_PATH} \
                --BQSR_PATH ${BQSR_PATH} --BQSR_RECAL_PATH ${BQSR_RECAL_PATH}  --FINAL_BAM_PATH ${FINAL_BAM_PATH} \
                --REF ${REF} --dbsnp ${dbsnp} --TMP_PATH ${TMP_PATH}")
            POSTBWA_JOBID=$(echo ${POSTBWA_JOB} | awk '{print $NF}')
            
            # Depth of coverage: PRE_BAM_PATH (RG가 안붙어서 안됨) vs AddOrReplaceReadGroups_PATH
            DOC_JOB=$(sbatch --dependency=afterok:${POSTBWA_JOBID} --job-name='doc_'${Sample_ID}"_"${TISSUE}"_"${hg} --output=${logPath}"/03.depthofcoverage/"${Sample_ID}"_"${TISSUE}".out" --error=${logPath}"/03.depthofcoverage/"${Sample_ID}"_"${TISSUE}".err" \
                --cpus-per-task=5 --mem-per-cpu=14GB --time=6:00:00 --partition=cpu --qos nstumor \
                --wrap="bash depthofcoverage.sh --REF ${REF}  --DOC_PATH ${DOC_PATH} --BAM_PATH ${AddOrReplaceReadGroups_PATH}  --INTERVAL ${INTERVAL}")
            DOC_JOBID=$(echo ${DOC_JOB} | awk '{print $NF}')

            # Mosdepth for PRE_BAM_PATH
            MOSDEPTH_JOB=$(sbatch --dependency=afterok:${POSTBWA_JOBID} --job-name='mosdepth_'${Sample_ID}"_"${TISSUE}"_"${hg} --output=${logPath}"/04.mosdepth/"${Sample_ID}"_"${TISSUE}".out" --error=${logPath}"/04.mosdepth/"${Sample_ID}"_"${TISSUE}".err" \
                --cpus-per-task=4 --mem-per-cpu=14GB --time=5:00:00 --partition=cpu --qos nstumor \
                --wrap="bash mosdepth.sh --THREADS 8 --INPUT_BAM ${PRE_BAM_PATH} --INTERVAL ${INTERVAL} --BAM_MOSDEPTH_PREFIX ${MOSDEPTH_DIR}"/"${Sample_ID}"_"${TISSUE}")
            MOSDEPTH_JOBID=$(echo ${MOSDEPTH_JOB} | awk '{print $NF}')

        fi    
    done
done


# qsub -pe smp 5 -e ${logPath}"/01.bwa" -o ${logPath}"/01.bwa" -N "bwa_"${Sample_ID}"_"${TISSUE}"_"${hg} -hold_jid 'FP_'${Sample_ID}"_"${TISSUE}",FQC_"${Sample_ID}"_"${TISSUE} bwamem_pipe_01.bwa.sh \
# --FASTQ_PATH_1 ${FASTQ_PATH_1} --FASTQ_PATH_2 ${FASTQ_PATH_2}  --PRE_BAM_PATH ${PRE_BAM_PATH} --REF ${REF} --TMP_PATH ${TMP_PATH}

# qsub -pe smp 8 -e ${logPath}"/02.postbwa" -o ${logPath}"/02.postbwa" -N "postbwa_"${Sample_ID}"_"${TISSUE}"_"${hg} -hold_jid 'bwa_'${Sample_ID}"_"${TISSUE}"_"${hg} bwamem_pipe_02.postbwa.sh \
# --Sample_ID ${Sample_ID}  --TISSUE ${TISSUE} --PRE_BAM_PATH ${PRE_BAM_PATH} --MarkDuplicate_PATH ${MarkDuplicate_PATH} --AddOrReplaceReadGroups_PATH ${AddOrReplaceReadGroups_PATH} \
# --BQSR_PATH ${BQSR_PATH} --BQSR_RECAL_PATH ${BQSR_RECAL_PATH}  --FINAL_BAM_PATH ${FINAL_BAM_PATH} \
# --REF ${REF} --dbsnp ${dbsnp} --TMP_PATH ${TMP_PATH}

# qsub -pe smp 5 -e ${logPath}"/03.depthofcoverage" -o ${logPath}"/03.depthofcoverage" -N "doc_"${Sample_ID}"_"${TISSUE}"_"${hg} -hold_jid 'bwa_'${Sample_ID}"_"${TISSUE}"_"${hg}',postbwa_'${Sample_ID}"_"${TISSUE}"_"${hg} depthofcoverage.sh \
# --REF ${REF}  --DOC_PATH ${DOC_PATH} --BAM_PATH ${AddOrReplaceReadGroups_PATH}  --INTERVAL ${INTERVAL}

# qsub -pe smp 4 -e $logPath"/04.mosdepth" -o $logPath"/04.mosdepth" -N "mosdepth_"${Sample_ID}"_"${TISSUE}"_"${hg} -hold_jid 'bwa_'${Sample_ID}"_"${TISSUE}"_"${hg}',postbwa_'${Sample_ID}"_"${TISSUE}"_"${hg}  "mosdepth.sh" \
#     --THREADS 8 \
#     --INPUT_BAM ${PRE_BAM_PATH} \
#     --INTERVAL ${INTERVAL} \
#     --BAM_MOSDEPTH_PREFIX ${MOSDEPTH_DIR}"/"${Sample_ID}"_"${TISSUE}
