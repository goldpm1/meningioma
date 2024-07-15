#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long REF:,INTERVAL:,INPUT_BAM_PATH:,OUTPUT_MPILEUP_PATH:, -- "$@")
then
    echo "ERROR: invalid options"
    exit 1
fi

eval set -- $options

while true; do
    case "$1" in
        -h|--help)
            echo "Usage"
        shift ;;
        --REF)
            REF=$2
        shift 2 ;;
        --INTERVAL)
            INTERVAL=$2
        shift 2 ;;
        --INPUT_BAM_PATH)
            INPUT_BAM_PATH=$2
        shift 2 ;;
        --OUTPUT_MPILEUP_PATH)
            OUTPUT_MPILEUP_PATH=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

#samtools mpileup -s -Q 20 -f ${REF} -l ${INTERVAL} -o ${OUTPUT_MPILEUP_PATH} ${INPUT_BAM_PATH}

echo -e "/opt/Yonsei/bcftools/1.17/bcftools-1.17/bcftools mpileup --indels-2.0 -Oz -a DP,AD \
    -R ${INTERVAL} -f ${REF}  -o ${OUTPUT_MPILEUP_PATH} ${INPUT_BAM_PATH}"

/opt/Yonsei/bcftools/1.17/bcftools-1.17/bcftools mpileup --indels-2.0 -Oz -a DP,AD \
    -R ${INTERVAL} -f ${REF}  -o ${OUTPUT_MPILEUP_PATH} ${INPUT_BAM_PATH}


#/opt/Yonsei/bcftools/1.17/bcftools-1.17/bcftools mpileup --indels-2.0 -Oz -a DP,AD     -R /home/goldpm1/resources/NF2.exome.proteincoding.bed -f /home/goldpm1/reference/genome.fa  -o /data/project/Meningioma/00.Targetted/07.mpileup/hg38/Dura/01.mpileup/230419_Dura.vcf.gz /data/project/Meningioma/02.Align/hg38/Dura/05.Final_bam/230419_Dura.bam