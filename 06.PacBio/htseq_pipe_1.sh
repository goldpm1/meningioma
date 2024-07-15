#!/bin/bash
#$ -cwd
#$ -S /bin/bash

INPUT=$1
OUTPUT=$2
GTFPath=$3

htseq-count -f bam -r name -t exon -i gene_id -a 20 -m intersection-strict --stranded=no  ${INPUT} ${GTFPath} > ${OUTPUT}




# # PATH #
# BAM_PATH=/home/jyhong906/Project/CRC/RNA/03.BAM
# HTseq_PATH=/home/jyhong906/Project/CRC/RNA/04.HTseq
# ​
# # REFERENCE #
# gtf=/data/resource/annotation/human/UCSC/hg38/Genes/genes.gtf
# ​
# # SAMPLE #
# library_name=$1
# ​
# htseq-count \
# -f bam \
# -r name \
# -s reverse \
# -a 10 \
# -t exon \
# -i gene_id \
# -m intersection-nonempty \
# $BAM_PATH/$library_name'_Aligned.out.bam' \
# $gtf > $HTseq_PATH/$library_name'.htseq.count.txt