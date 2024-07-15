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

for sublog in 52.isoseq_refine  53.isoseq_cluster 54.pbmm2 55.isoseq_collapse 56.pigeon; do  #01.WGRS
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done





# 2. Isoseq-RNA
FASTQ_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/00.raw"
UBAM_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/00.raw/ubam"
UBAM_REFINE_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/01.QC/01.refine"
UBAM_CLUSTER_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/01.QC/02.cluster"
BAM_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/02.Align"
COLLAPSE_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/02.Align/collapse"
PIGEON_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/03.pigeon"
MATRIX_DIR=${PROJECT_DIR}"/01.PacBio/Isoseq_RNA/03.matrix"
PRIMER_XML=${PROJECT_DIR}"/01.PacBio/Iso-Seq_v2_Barcoded_cDNA_Primers.barcodeset.xml"
PRIMER_FASTA=${PROJECT_DIR}"/01.PacBio/IsoSeq_v2_primers_12.fasta"

TMP_PATH=${UBAM_DIR}"/temp"
for subdir in ${UBAM_REFINE_DIR} ${UBAM_CLUSTER_DIR} ${BAM_DIR} ${COLLAPSE_DIR} ${MATRIX_DIR}; do 
    if [ ! -d $subdir ] ; then
        mkdir -p $subdir
    fi
done



FASTQ=$(find "$FASTQ_DIR" -type f -name "*.fastq.gz")
FASTQ_LIST=(${FASTQ// / })                   # 이를 배열 (list)로 만듬
UBAM=$(find "$UBAM_DIR" -type f -name "*.bam" | grep -v 'pbi')
UBAM_LIST=(${UBAM// / })                   # 이를 배열 (list)로 만듬

PRESET="ISOSEQ"

for idx in ${!UBAM_LIST[@]}              # @ 배열의 모든 element    #! : indexing
do
        UBAM_PATH=${UBAM_LIST[idx]}         # idx번째의 파일명을 담아둔다
        ID=${UBAM_LIST[idx]/"${UBAM_DIR}/"/}  
        ID=${ID%".bam"}
        #echo -e "FASTQ_PATH: "$FASTQ_PATH   #  /data/project/Meningioma/01.PacBio/Isoseq_RNA/00.raw/m84213_240222_061040_s2.IsoSeqX_bc03_5p--IsoSeqX_3p.hifi_reads.fastq.gz
        echo -e "ID: "$ID   #  m84213_240222_061040_s2.IsoSeqX_bc03_5p--IsoSeqX_3p.hifi_reads 


        # [52.Isoseq_refine:   generate FLNC bam ]
        #OUTPUT_BAM=${UBAM_REFINE_DIR}"/"${ID}"_refine.bam
        # qsub -pe smp 2 -e $logPath"/52.isoseq_refine" -o $logPath"/52.isoseq_refine" -N 'PacBio_52.isoseq_refine.'${ID}  ${CURRENT_PATH}"/PacBio_IsoSeq_pipe_52.isoseq_refine.sh" \
        #     --INPUT_BAM ${UBAM_PATH} --PRIMER_XML ${PRIMER_XML} --OUTPUT_BAM ${OUTPUT_BAM}

        # [53.Isoseq_cluster: generate unpolished(clustered) bam ]
        # FILE_LIST=${UBAM_CLUSTER_DIR}"/flnc.fofn"
        # INPUT_BAM=${UBAM_REFINE_DIR}"/"${ID}"_refine.bam"
        # OUTPUT_BAM=${UBAM_CLUSTER_DIR}"/"${ID}"_clusterd.bam"
        # qsub -pe smp 10 -e $logPath"/53.isoseq_cluster" -o $logPath"/53.isoseq_cluster" -N 'PacBio_53.isoseq_cluster' ${CURRENT_PATH}"/PacBio_IsoSeq_pipe_53.isoseq_cluster.sh" \
        #     --FILE_LIST ${FILE_LIST}  --INPUT_BAM ${INPUT_BAM} --OUTPUT_BAM ${OUTPUT_BAM}

        # [43. pbmm : Aligner]
        # qsub -pe smp 2 -e $logPath"/02.Isoseq" -o $logPath"/02.Isoseq" -N 'PacBio_01.'${ID}  ${CURRENT_PATH}"/PacBio_pipe_01.pbmm2_mkindex.sh" \
        #     --REF ${REF} --REF_MMI "/home/goldpm1/reference/genome."${PRESET}".mmi" --PRESET "ISOSEQ"
        # INPUT_BAM=${UBAM_CLUSTER_DIR}"/"${ID}"_clusterd.bam"
        # qsub -pe smp 8 -e $logPath"/54.pbmm2" -o $logPath"/54.pbmm2" -N 'PacBio_54.pbmm2.'${ID}  ${CURRENT_PATH}"/PacBio_pipe_01.pbmm2.sh" \
        #     --REF_MMI "/home/goldpm1/reference/genome.CCS/genome.CCS.mmi" \
        #     --BAM_PATH ${INPUT_BAM} \
        #     --PRESET "ISOSEQ" \
        #     --OUTPUT_BAM ${BAM_DIR}"/"${ID}".52.bam" \
        #     --SAMPLENAME ${ID}

        # [55. Isoseq_collapse]
        # INPUT_BAM=${BAM_DIR}"/"${ID}".52.bam"
        # OUTPUT_GFF=${COLLAPSE_DIR}"/"${ID}".gff"
        # qsub -pe smp 5 -e $logPath"/55.isoseq_collapse" -o $logPath"/55.isoseq_collapse" -N 'PacBio_55.isoseq_collapse.'${ID}  ${CURRENT_PATH}"/PacBio_IsoSeq_pipe_55.isoseq_collapse.sh" \
        #     --INPUT_BAM ${INPUT_BAM} \
        #     --OUTPUT_GFF ${OUTPUT_GFF}

        # [56. pigeon]
        GENCODE_GTF="/home/goldpm1/resources/gencode.v38.annotation.gtf"
        #GENCODE_ANNOTATION_GTF="/home/goldpm1/resources/gencode.v38.annotation.gtf"
        INPUT_GFF=${COLLAPSE_DIR}"/"${ID}".gff"
        OUTPUT_DIR=${PIGEON_DIR}
        qsub -pe smp 5 -e $logPath"/56.pigeon" -o $logPath"/56.pigeon" -N 'PacBio_56.pigeon.'${ID}  ${CURRENT_PATH}"/PacBio_IsoSeq_pipe_56.pigeon.sh" \
            --GENCODE_GTF ${GENCODE_GTF} \
            --REF ${REF} \
            --ID ${ID} \
            --INPUT_GFF ${INPUT_GFF} \
            --OUTPUT_DIR ${OUTPUT_DIR}


        # HTseq : 이렇게 하면 곤란하다
        #GTFPath="/home/goldpm1/resources/gencode.v38.primary_assembly.annotation.gtf"     # Ensembl ID로 나옴
        # GTFPath="/data/resource/annotation/human/UCSC/hg38/Genes/genes.gtf "       # Gene symbol로 나옴
        # qsub -pe smp 5 -o ${logPath}"/02.Isoseq" -e ${logPath}"/02.Isoseq" -N "HT_"${ID}  -hold_jid 'PacBio_01.'${ID} ${CURRENT_PATH}"/htseq_pipe_1.sh" \
        #     ${BAM_DIR}"/"${ID}".bam" \
        #     ${MATRIX_DIR}"/"${ID}".tsv" \
        #     ${GTFPath}

done


# ls ${UBAM_REFINE_DIR}"/*_Tumor_refine .bam"  > ${UBAM_CLUSTER_DIR}"/flnc.fofn"



# [53.Isoseq_clustering ]
