#!/bin/bash
#$ -cwd
#$ -S /bin/bash

idx=$1
AnalysisPath=$2
MT_Path=$3
MFSNP=$4
MFIND=$5

Input_mutect=${MT_Path}"/CR."$idx".mutect2.filtered.vcf"
Input_mutect_exon=${MT_Path}"/CR."$idx".mutect2.filtered.exon.vcf"


# SNP

python3 Mosaic_pipe_4.formatchange.py $MFSNP
Input_MF_SNP=$MFSNP".bed"

Output_MF_SNP=$AnalysisPath"/CR."$idx"/CR."$idx".predictions.filtered.snp.final.vcf"
Output_MF_SNP_exon=$AnalysisPath"/CR."$idx"/CR."$idx".predictions.filtered.exon.snp.final.vcf"

bedtools intersect -header -a $Input_mutect -b $Input_MF_SNP > $Output_MF_SNP
bedtools intersect -header -a $Input_mutect_exon -b $Input_MF_SNP > $Output_MF_SNP_exon

echo "SNP done"


# IND
python3 Mosaic_pipe_4.formatchange.py $MFIND
Input_MF_IND=$MFIND".bed"

Output_MF_IND=$AnalysisPath"/ind.CR."$idx"/CR."$idx".predictions.filtered.indel.final.vcf"
Output_MF_IND_exon=$AnalysisPath"/ind.CR."$idx"/CR."$idx".predictions.filtered.exon.indel.final.vcf"

bedtools intersect -header -a $Input_mutect -b $Input_MF_IND > $Output_MF_IND
bedtools intersect -header -a $Input_mutect_exon -b $Input_MF_IND > $Output_MF_IND_exon

echo "IND done"


# 합치기
Output_MF=$AnalysisPath"/CR."$idx"/CR."$idx".predictions.filtered.final.vcf"
Output_MF_exon=$AnalysisPath"/CR."$idx"/CR."$idx".predictions.filtered.exon.final.vcf"

cat $Output_MF_SNP > $Output_MF || grep -v '#' $Output_MF_IND >> $Output_MF
cat $Output_MF_SNP_exon > $Output_MF_exon || grep -v '#' $Output_MF_IND_exon >> $Output_MF_exon

echo "Merge done"



gzip -f $Output_MF_SNP
gzip -f $Output_MF_SNP_exon

gzip -f $Output_MF_IND
gzip -f $Output_MF_IND_exon

gzip -f $Output_MF
gzip -f $Output_MF_exon

echo "gzip done"
