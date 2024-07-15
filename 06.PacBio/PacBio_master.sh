#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

PROJECT_DIR="/data/project/Meningioma"

PON="/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
REF="/home/goldpm1/reference/genome.fa"
hg="hg38"
gnomad="/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi

for sublog in 01.WGRS 02.Isoseq; do  #01.WGRS
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done



BAM_DIR=${PROJECT_DIR}"/01.PacBio/WGRS/02.Align"
TMP_PATH=${BAM_DIR}"/temp"
VCF_DIR=${PROJECT_DIR}"/01.PacBio/WGRS/03.vcf"
SV_DIR=${PROJECT_DIR}"/01.PacBio/WGRS/04.sv"
for subdir in ${BAM_DIR} ${VCF_DIR} ${SV_DIR} ${TMP_PATH}; do 
    if [ ! -d $subdir ] ; then
        mkdir -p $subdir
    fi
done

BAM=$(find "$BAM_DIR" -type f -name "*.bam")
BAM_LIST=(${BAM// / })                   # 이를 배열 (list)로 만듬

PRESET="CCS"  # "ISOSEQ"

for idx in ${!BAM_LIST[@]}              # @ 배열의 모든 element    #! : indexing
do
        BAM_PATH=${BAM_LIST[idx]}         # idx번째의 파일명을 담아둔다
        ID=${BAM_LIST[idx]/"${BAM_DIR}/"/}  
        ID=${ID%".bam"}
        echo -e "BAM_PATH: "$BAM_PATH   # /data/project/Meningioma/01.PacBio/WGRS/00.raw/m84213_240222_041127_s1.hifi_reads.bc2010.bam
        echo -e "ID: "$ID   # m84213_240222_041127_s1.hifi_reads.bc2010


        # [pbmm : Aligner]   한번만 하면 됨
        # qsub -pe smp 2 -e $logPath"/01.WGRS" -o $logPath"/01.WGRS" -N 'PacBio_01.'${ID}  ${CURRENT_PATH}"/PacBio_pipe_01.pbmm2_mkindex.sh" \
        #     --REF ${REF} --REF_MMI "/home/goldpm1/reference/genome."${PRESET}".mmi" --PRESET ${PRESET}

        # qsub -pe smp 8 -e $logPath"/01.WGRS" -o $logPath"/01.WGRS" -N 'PacBio_01.'${ID}  ${CURRENT_PATH}"/PacBio_pipe_01.pbmm2.sh" \
        #     --REF_MMI "/home/goldpm1/reference/genome."${PRESET}".mmi" \
        #     --BAM_PATH ${BAM_PATH} \
        #     --PRESET ${PRESET} \
        #     --OUTPUT_BAM "/data/project/Meningioma/01.PacBio/WGRS/02.Align/230526_Tumor_Pacbio.bam" \
        #     --SAMPLENAME "230526_Tumor"


        # [variant caller]
        OUTPUT_VCF_PATH="/data/project/Meningioma/01.PacBio/WGRS/03.vcf/230526_Tumor_Pacbio.vcf"
        OUTPUT_MUTECT_GZ="/data/project/Meningioma/01.PacBio/WGRS/03.vcf/230526_Tumor_Pacbio.mutect2.vcf.gz"
        OUTPUT_MUTECT="/data/project/Meningioma/01.PacBio/WGRS/03.vcf/230526_Tumor_Pacbio.mutect2.vcf"
        OUTPUT_FMC_PATH="/data/project/Meningioma/01.PacBio/WGRS/03.vcf/230526_Tumor_Pacbio.mutect2.FMC.vcf"
        OUTPUT_FMC_HF_PATH="/data/project/Meningioma/01.PacBio/WGRS/03.vcf/230526_Tumor_Pacbio.mutect2.FMC.HF.vcf"
        OUTPUT_FMC_HF_RMBLACK_PATH="/data/project/Meningioma/01.PacBio/WGRS/03.vcf/230526_Tumor_Pacbio.mutect2.FMC.HF.RMBLACK.vcf"
        #"/home/goldpm1/miniconda3/envs/PacBio/lib/python3.7/site-packages/pbcore/data/lambdaNEB.fa" 
        SAMPLE_THRESHOLD="all"   # "all"
        DP_THRESHOLD=10
        ALT_THRESHOLD=1
        REMOVE_MULTIALLELIC="True"
        PASS="True"
        REMOVE_MITOCHONDRIAL_DNA="True"
        BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
        # qsub -pe smp  5 -e $logPath"/01.WGRS" -o $logPath"/01.WGRS" -N 'PacBio_02.'${ID}  -hold_jid 'PacBio_01.'${ID}  ${CURRENT_PATH}"/PacBio_pipe_02.gcpp.sh" \
        #     --BAM_PATH "/data/project/Meningioma/01.PacBio/WGRS/02.Align/230526_Tumor_Pacbio.bam" \
        #     --REF_MMI "/home/goldpm1/reference/genome."${PRESET}".mmi" --REF ${REF} --PON ${PON} --gnomad ${gnomad} --TMP_PATH ${TMP_PATH} \
        #     --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} --BLACKLIST ${BLACKLIST}  \
        #     --OUTPUT_VCF_PATH ${OUTPUT_VCF_PATH} --OUTPUT_MUTECT_GZ ${OUTPUT_MUTECT_GZ} --OUTPUT_MUTECT ${OUTPUT_MUTECT} \
        #     --OUTPUT_FMC_PATH ${OUTPUT_FMC_PATH} --OUTPUT_FMC_HF_PATH ${OUTPUT_FMC_HF_PATH} --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH}

        # # [SV caller]
        # qsub -pe smp 8 -e $logPath"/01.WGRS" -o $logPath"/01.WGRS" -N 'PacBio_03.'${ID}  -hold_jid 'PacBio_01.'${ID} ${CURRENT_PATH}"/PacBio_pipe_03.pbsv.sh" \
        #     --BAM_PATH ${BAM_PATH} \
        #     --OUTPUT_SIG ${SV_DIR}"/230526_Tumor_Pacbio.svsig.gz" \
        #     --REF ${REF} \
        #     --OUTPUT_VCF ${SV_DIR}"/230526_Tumor_Pacbio.sv.vcf"


        # [Methylation calling:]   Jasmine은 해줬으니 pb-CpG-tools만 돌리면 된다
        qsub -pe smp  8 -e $logPath"/01.WGRS" -o $logPath"/01.WGRS" -N 'PacBio_04.'${ID}   ${CURRENT_PATH}"/PacBio_pipe_04.pb_CpG.sh" \
            --BAM_PATH "/data/project/Meningioma/01.PacBio/WGRS/02.Align/230526_Tumor_Pacbio.bam" \
            --OUTPUT_PREFIX "230526_Tumor_Pacbio_pb_CpG" \
            --MODEL "/home/goldpm1/tools/PacBio/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite" \
            --THREADS 16
            


done





# 2. Isoseq-RNA
# FASTQ_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/00.raw"
# TMP_PATH=${FASTQ_DIR}"/temp"
# BAM_DIR="/data/project/Meningioma/01.PacBio/Isoseq_RNA/02.Align"
# MATRIX_DIR="/data/project/Meningioma/01.PacBio/Isoseq_RNA/03.matrix"
# for subdir in ${BAM_DIR} ${MATRIX_DIR}; do 
#     if [ ! -d $subdir ] ; then
#         mkdir -p $subdir
#     fi
# done



# FASTQ=$(find "$FASTQ_DIR" -type f -name "*.fastq.gz")
# FASTQ_LIST=(${FASTQ// / })                   # 이를 배열 (list)로 만듬

# PRESET="ISOSEQ"

# for idx in ${!FASTQ_LIST[@]}              # @ 배열의 모든 element    #! : indexing
# do
#         FASTQ_PATH=${FASTQ_LIST[idx]}         # idx번째의 파일명을 담아둔다
#         ID=${FASTQ_LIST[idx]/"${FASTQ_DIR}/"/}  
#         ID=${ID%".hifi_reads.fastq.gz"}
#         #echo -e "FASTQ_PATH: "$FASTQ_PATH   #  /data/project/Meningioma/01.PacBio/Isoseq_RNA/00.raw/m84213_240222_061040_s2.IsoSeqX_bc03_5p--IsoSeqX_3p.hifi_reads.fastq.gz
#         echo -e "ID: "$ID   #  m84213_240222_061040_s2.IsoSeqX_bc03_5p--IsoSeqX_3p.hifi_reads 


#         # [pbmm : Aligner]
#         # qsub -pe smp 2 -e $logPath"/02.Isoseq" -o $logPath"/02.Isoseq" -N 'PacBio_01.'${ID}  ${CURRENT_PATH}"/PacBio_pipe_01.pbmm2_mkindex.sh" \
#         #     --REF ${REF} --REF_MMI "/home/goldpm1/reference/genome."${PRESET}".mmi" --PRESET ${PRESET}

#         # qsub -pe smp 8 -e $logPath"/02.Isoseq" -o $logPath"/02.Isoseq" -N 'PacBio_01.'${ID}  ${CURRENT_PATH}"/PacBio_pipe_01.pbmm2.sh" \
#         #     --REF_MMI "/home/goldpm1/reference/genome."${PRESET}".mmi" \
#         #     --BAM_PATH ${FASTQ_PATH} \
#         #     --PRESET ${PRESET} \
#         #     --OUTPUT_BAM ${BAM_DIR}"/"${ID}".bam" \
#         #     --SAMPLENAME ${ID}

#         #GTFPath="/home/goldpm1/resources/gencode.v38.primary_assembly.annotation.gtf"     # Ensembl ID로 나옴
#         GTFPath="/data/resource/annotation/human/UCSC/hg38/Genes/genes.gtf "       # Gene symbol로 나옴
#         qsub -pe smp 5 -o ${logPath}"/02.Isoseq" -e ${logPath}"/02.Isoseq" -N "HT_"${ID}  -hold_jid 'PacBio_01.'${ID} ${CURRENT_PATH}"/htseq_pipe_1.sh" \
#             ${BAM_DIR}"/"${ID}".bam" \
#             ${MATRIX_DIR}"/"${ID}".tsv" \
#             ${GTFPath}

# done
