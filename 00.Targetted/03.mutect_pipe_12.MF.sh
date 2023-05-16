#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long REF:,BAM_DIR:,CASE_BAM_PATH:,ID:,MT_SNP_VCF:,MT_INDEL_VCF:,MF_DIR:,F_bed:,F_bed_ind:,MF_TOOL_DIR:, -- "$@")
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
        --BAM_DIR)
            BAM_DIR=$2
        shift 2 ;;
        --CASE_BAM_PATH)
            CASE_BAM_PATH=$2
        shift 2 ;;
        --ID)
            ID=$2
        shift 2 ;;
        --MT_SNP_VCF)
            MT_SNP_VCF=$2
        shift 2 ;;
        --MT_INDEL_VCF)
            MT_INDEL_VCF=$2
        shift 2 ;;
        --MF_DIR)
            MF_DIR=$2
        shift 2 ;;
        --F_bed)
            F_bed=$2
        shift 2 ;;
        --F_bed_ind)
            F_bed_ind=$2
        shift 2 ;;
        --MF_TOOL_DIR)
            MF_TOOL_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


# REF="/home/goldpm1/reference/genome.fa"
# BAM_DIR="/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/05.Final_bam"
# F_bed="/data/project/MRS/Resource/MF/fout_snv.bed"
# F_bed_ind="/data/project/MRS/Resource/MF/fout_ind.bed"
# MF_TOOL_DIR="/home/goldpm1/tools/MosaicForecast2"



MF_DIR_01=${MF_DIR}"/01.input"
MF_DIR_02=${MF_DIR}"/02.output"

for folder in  ${MF_DIR_01}  ${MF_DIR_02}  ; do
    if [ ! -d $folder ] ; then
        mkdir $folder
    fi
done



######################### [SNP] #########################

  # 01. MT_SNP_VCF의 특정 구간을 filter out
  bedtools intersect -v -a ${MT_SNP_VCF} -b ${F_bed} > ${MF_DIR_01}"/"$ID".MT.filtered.snv.vcf"

  # 02. VCF -> BED
  python3 ${MF_TOOL_DIR}"/argparse_YS/MFscript_convert_VCFtoBED.mx.py" \
    --INPUT_VCF  ${MF_DIR_01}"/"$ID".MT.filtered.snv.vcf" \
    --OUTPUT_BED ${MF_DIR_01}"/"$ID".MT.filtered.snv.bed" \
    --ID ${ID}

  # 03.  참고 : $BAM_DIR"/"$ID".bam" 형식으로 존재해야 함
  #mkdir ${MF_DIR_02}"/"$ID
  python3 ${MF_TOOL_DIR}"/argparse_YS/ReadLevel_Features_extraction.py" \
        --INPUT_BED  ${MF_DIR_01}"/"$ID".MF.filtered.snv.bed" \
        --OUTPUT        ${MF_DIR_02}"/"$ID".features" \
        --BAM_DIR       ${BAM_DIR} \
        --REF                 ${REF} \
        --UNIMAP_MAPPABILITY_BIGWIGFILE "/data/project/MRS/Resource/MF/k24.umap.wg.bw" \
        --NUMBERJOBS 2 \
        --SEQ_FILE_FORMAT "bam"

#   python3 "/home/goldpm1/tools/MosaicForecast2/argparse_YS/ReadLevel_Features_extraction.py" \
#         --INPUT_BED  "/data/project/Meningioma/00.Targetted/06.Mutect/12.MF/01.input/230106_Dura.MF.filtered.snv.bed" \
#         --OUTPUT        "/data/project/Meningioma/00.Targetted/06.Mutect/12.MF/02.output/230106_Dura.features" \
#         --BAM_DIR       "02.Align/hg38/Dura/05.Final_bam" \
#         --REF                 "/home/goldpm1/reference/genome.fa" \
#         --UNIMAP_MAPPABILITY_BIGWIGFILE "/data/project/MRS/Resource/MF/k24.umap.wg.bw" \
#         --NUMBERJOBS 2 \
#         --SEQ_FILE_FORMAT "bam"


# 04.
#   python3 ${MF_TOOL_DIR}"/argpare_YS/MFmodi.py" ${MF_DIR_02}"/"$ID"/"$ID


#   05. Prediction SNP
#   Rscript  ${MF_TOOL_DIR}"/Prediction.R" \
#             ${MF_DIR_02}"/"$ID"/"$ID".features.modi" \
#             ${MF_TOOL_DIR}"/models_trained/50xRFmodel_addRMSK_Refine.rds" \
#             Refine \
#             ${MF_DIR_02}"/"$ID"/"$ID".SNP.Predictions"
#   ###
#   cat ${MF_DIR_02}"/"$ID"/"$ID".SNP.Predictions" | grep mosaic | grep SNP > ${MF_DIR_02}"/"$ID"/"$ID".SNP.Predictions.mosaic"

#   echo "SNP done"
#   date



# # [INDEL]

#   # Filter 한 vcf,  그걸 가공한 bed 파일을 {MF_DIR_01}에 넣는다

#   bedtools intersect -v -a $MT_Path"/"$ID".mutect2.filtered.exon.indel.vcf" -b $F_bed_ind > ${MF_DIR_01}"/"$ID".MF.Filterout.indel.vcf"

#   python $MF_TOOL_DIR"/MFscript_convert_VCFtoBED.mx.py" $ID ${MF_DIR_01}  ${MF_DIR_01} 'indel'

#   ###
#   mkdir ${MF_DIR_02}"/ind."$ID
#   python3 $MF_TOOL_DIR"/ReadLevel_Features_extraction.py" \
#           ${MF_DIR_01}"/"$ID".MF.filtered.indel.input" \
#           ${MF_DIR_02}"/ind."$ID"/"$ID".features" \
#           $BAM_DIR \
#           "/data/project/MRS/Resource/reference/genome.fa"\
#           "/data/project/MRS/Resource/MF/k24.umap.wg.bw"  2 bam

#   python /data/project/MRS/script/MFmodi.py ${MF_DIR_02}"/ind."$ID"/"$ID


#   ########Prediction INS
#   Rscript  $MF_TOOL_DIR"/Prediction.R" \
#   ${MF_DIR_02}"/ind."$ID"/"$ID".features.modi" \
#   $MF_TOOL_DIR"/models_trained/deletions_250x.RF.rds" \
#   Refine \
#   ${MF_DIR_02}"/ind."$ID"/"$ID".DEL.Predictions"
#   cat ${MF_DIR_02}"/ind."$ID"/"$ID".DEL.Predictions" | egrep 'hap=3' | grep DEL > ${MF_DIR_02}"/ind."$ID"/"$ID".DEL.Predictions.mosaic"

#   echo "INS done"

#   ########Prediction DEL
#   Rscript  $MF_TOOL_DIR"/Prediction.R" \
#   ${MF_DIR_02}"/ind."$ID"/"$ID".features.modi" \
#   $MF_TOOL_DIR"/models_trained/insertions_250x.RF.rds" \
#   Refine \
#   ${MF_DIR_02}"/ind."$ID"/"$ID".INS.Predictions"
#   cat ${MF_DIR_02}/"ind."$ID"/"$ID".INS.Predictions" | egrep 'hap=3' | grep INS > ${MF_DIR_02}"/ind."$ID"/"$ID".INS.Predictions.mosaic"

#   echo "DEL done"

# cat ${MF_DIR_02}"/ind."$ID"/"$ID".DEL.Predictions.mosaic" ${MF_DIR_02}"/ind."$ID"/"$ID".INS.Predictions.mosaic" > ${MF_DIR_02}"/ind."$ID"/"$ID".IND.Predictions.mosaic"
