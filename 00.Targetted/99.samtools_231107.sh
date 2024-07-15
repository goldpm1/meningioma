#!/bin/bash
#$ -cwd
#$ -S /bin/bash

#PATH#
TMP_PATH="/data/project/Meningioma/00.Targetted/temp"

#REF="/data/resource/reference/human/NCBI/GRCh38_GATK/BWAIndex/genome.fa"
REF="/home/goldpm1/reference/genome.fa"
#INTERVAL=/data/project/MultiSampleMosaicCaller/MRS_dataset/resource/Interval/for_samtools_mpileup.bed


for date in 230822; do
    echo ${date}
    INTERVAL="/data/project/Meningioma/00.Targetted/07.mpileup/hg38/Dura/00.interval/"${date}".interval.bed"
    # WES
    INPUT_PATH="/data/project/Meningioma/02.Align/hg38/Dura/05.Final_bam/"${date}"_Dura.bam"
    OUTPUT_PATH="/data/project/Meningioma/00.Targetted/07.mpileup/hg38/Dura/01.mpileup/"${date}"_WES.pileup"
    echo "1.WES : /opt/Yonsei/samtools/1.7/samtools mpileup -f ${REF} -l ${INTERVAL} -o ${OUTPUT_PATH} ${INPUT_PATH}"
    /opt/Yonsei/samtools/1.7/samtools mpileup -f ${REF} -l ${INTERVAL} -o ${OUTPUT_PATH} ${INPUT_PATH}
    echo -e "\n"

    # Amplicon (original)
    INPUT_PATH="/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/01.Pre_bam/"${date}".sorted.bam"
    OUTPUT_PATH="/data/project/Meningioma/00.Targetted/07.mpileup/hg38/Dura/01.mpileup/"${date}"_amplicon_duplicated.pileup"
    echo "2. Amplicon (duplicated) : /opt/Yonsei/samtools/1.7/samtools mpileup -f ${REF} -l ${INTERVAL} -o ${OUTPUT_PATH} ${INPUT_PATH}"
    /opt/Yonsei/samtools/1.7/samtools mpileup -f ${REF} -l ${INTERVAL} -o ${OUTPUT_PATH} ${INPUT_PATH}
    echo -e "\n"

    # Amplicon (deduplicated)
    INPUT_PATH="/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/05.Final_bam/"${date}"_Dura.bam"
    OUTPUT_PATH="/data/project/Meningioma/00.Targetted/07.mpileup/hg38/Dura/01.mpileup/"${date}"_amplicon_deduplicated.pileup"
    echo "3. Amplicon (deduplicated) :  /opt/Yonsei/samtools/1.7/samtools mpileup -f ${REF} -l ${INTERVAL} -o ${OUTPUT_PATH} ${INPUT_PATH}"
    /opt/Yonsei/samtools/1.7/samtools mpileup -f ${REF} -l ${INTERVAL} -o ${OUTPUT_PATH} ${INPUT_PATH}
    echo -e "\n\n"
done