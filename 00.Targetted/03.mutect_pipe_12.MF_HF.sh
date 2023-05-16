#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ID=$1
AnalysisPath=$2

# AnalysisPath="/data/project/craniosynostosis/9.mosaic/analysis_HF/"

REF="/home/goldpm1/reference/genome.fa"
F_bed="/data/project/MRS/Resource/MF/fout_snv.bed"
F_bed_ind="/data/project/MRS/Resource/MF/fout_ind.bed"

HC_Path="/data/project/craniosynostosis/ploidy/ploidy4"
InputPath="/data/project/craniosynostosis/9.mosaic/input_HF"
BamPath="/data/project/craniosynostosis/0.raw/0.bam"
#MF_Path="/opt/Yonsei/MosaicForecast/0.0.1"
MF_Path="/home/goldpm1/tools/MosaicForecast2"
echo $ID


# [SNP]

# Filter 한 vcf,  그걸 가공한 bed 파일을 InputPath에 넣는다

bedtools intersect -v -a $HC_Path"/"$ID".HF.snpindel.ploidy4.remdup.remblack.snp.vcf" -b $F_bed > $InputPath"/"$ID".MF.Filterout.snv.vcf"

python $MF_Path"/MFscript_convert_VCFtoBED.mx.py" $ID $InputPath $InputPath 'snv'

####
mkdir $AnalysisPath"/"$ID
python3 $MF_Path"/ReadLevel_Features_extraction.py" \
$InputPath"/"$ID".MF.filtered.snv.input" \
$AnalysisPath"/"$ID"/"$ID".features" \
$BamPath \
"/data/project/MRS//Resource/reference/genome.fa" \
"/data/project/MRS/Resource/MF/k24.umap.wg.bw"  2 bam
# ####
# python $MF_Path"/MFmodi.py" $AnalysisPath"/"$ID"/"$ID
# ######
# ########Prediction SNP
# Rscript  $MF_Path"/Prediction.R" \
# $AnalysisPath/$ID/$ID'.features.modi' \
# $MF_Path/models_trained/250xRFmodel_addRMSK_Refine.rds \
# Refine \
# $AnalysisPath/$ID/$ID'.SNP.Predictions'
# ###
# cat $AnalysisPath/$ID/$ID'.SNP.Predictions' | grep mosaic | grep SNP > $AnalysisPath/$ID/$ID'.SNP.Predictions.mosaic'
# ##
# ##
# ##


# [INDEL]
##
# bedtools intersect -v -a $HC_Path"/"$ID".HF.snpindel.ploidy4.remdup.remblack.indel.vcf" -b $F_bed_ind > $InputPath"/"$ID'.MF.Filterout.ind.vcf'
# ##
# python /data/project/MRS/script/MFscript_convert_VCFtoBED.mx.py $ID $InputPath   $InputPath 'ind'
##


mkdir $AnalysisPath'ind.'$ID/
python3 $MF_Path/ReadLevel_Features_extraction.py \
$InputPath/mx/$ID'.MF.filtered.ind.input' \
$AnalysisPath'ind.'$ID/$ID'.features' \
$DataPath \
$REF \
/data/project/MRS/Resource/MF/k24.umap.wg.bw  2 bam
##
# python /data/project/MRS/script/MFmodi.py $AnalysisPath/'ind.'$ID/$ID
# ##
# Rscript  $MF_Path/Prediction.R \
# $AnalysisPath/'ind.'$ID/$ID'.features.modi' \
# $MF_Path/models_trained/deletions_250x.RF.rds \
# Refine \
# $AnalysisPath/'ind.'$ID/$ID'.DEL.Predictions'
# cat $AnalysisPath/'ind.'$ID/$ID'.DEL.Predictions' | egrep 'hap=3' | grep DEL > $AnalysisPath/'ind.'$ID/$ID'.DEL.Predictions.mosaic'
# ######
# ######
# ######
# Rscript  $MF_Path/Prediction.R \
# $AnalysisPath/'ind.'$ID/$ID'.features.modi' \
# $MF_Path/models_trained/insertions_250x.RF.rds \
# Refine \
# $AnalysisPath/'ind.'$ID/$ID'.INS.Predictions'
# ######
# cat $AnalysisPath/'ind.'$ID/$ID'.INS.Predictions' | egrep 'hap=3' | grep INS > $AnalysisPath/'ind.'$ID/$ID'.INS.Predictions.mosaic'
#
# cat $AnalysisPath/'ind.'$ID/$ID'.DEL.Predictions.mosaic' $AnalysisPath/'ind.'$ID/$ID'.INS.Predictions.mosaic' > $AnalysisPath/'ind.'$ID/$ID'.IND.Predictions.mosaic'
