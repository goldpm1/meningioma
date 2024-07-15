#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if ! options=$(getopt -o h --long INPUT_VCF:,OUTPUT_VCF:,OUTPUT_BAMSNAP_DIR:,SAMPLE_ID:,BAM_DIR_LIST:,TITLE_LIST:,THRESHOLD_DEPTH:,WES_TUMOR_BED:,WES_DURA_VCF:,BCFTOOLS_MERGE_TXT:,BCFTOOLS_MERGE_OUTPUT_VCF_GZ:,OUTPUT_HEATMAP_PATH:, -- "$@")
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
        --INPUT_VCF)
            INPUT_VCF=$2
        shift 2 ;;
        --OUTPUT_VCF)
            OUTPUT_VCF=$2
        shift 2 ;;
        --OUTPUT_BAMSNAP_DIR)
            OUTPUT_BAMSNAP_DIR=$2
        shift 2 ;;
        --SAMPLE_ID)
            SAMPLE_ID=$2
        shift 2 ;;
        --BAM_DIR_LIST)
            BAM_DIR_LIST=$2
        shift 2 ;;
        --TITLE_LIST)
            TITLE_LIST=$2
        shift 2 ;;
        --THRESHOLD_DEPTH)
            THRESHOLD_DEPTH=$2
        shift 2 ;;
        --WES_TUMOR_BED)
            WES_TUMOR_BED=$2
        shift 2 ;;
        --WES_DURA_VCF)
            WES_DURA_VCF=$2
        shift 2 ;;
        --BCFTOOLS_MERGE_TXT)
            BCFTOOLS_MERGE_TXT=$2
        shift 2 ;;
        --BCFTOOLS_MERGE_OUTPUT_VCF_GZ)
            BCFTOOLS_MERGE_OUTPUT_VCF_GZ=$2
        shift 2 ;;
        --OUTPUT_HEATMAP_PATH)
            OUTPUT_HEATMAP_PATH=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo -e "python3 06.HC_pipe_05.pick_somatic.py --INPUT_VCF ${INPUT_VCF} --OUTPUT_VCF ${OUTPUT_VCF}  --OUTPUT_BAMSNAP_DIR ${OUTPUT_BAMSNAP_DIR} --SAMPLE_ID ${SAMPLE_ID} --BAM_DIR_LIST ${BAM_DIR_LIST} --TITLE_LIST ${TITLE_LIST} --THRESHOLD_DEPTH ${THRESHOLD_DEPTH} --WES_TUMOR_BED ${WES_TUMOR_BED} --WES_DURA_VCF ${WES_DURA_VCF}"
python3 "06.HC_pipe_05.pick_somatic.py" --INPUT_VCF ${INPUT_VCF} --OUTPUT_VCF ${OUTPUT_VCF}  --OUTPUT_BAMSNAP_DIR ${OUTPUT_BAMSNAP_DIR} --SAMPLE_ID ${SAMPLE_ID} --BAM_DIR_LIST ${BAM_DIR_LIST} --TITLE_LIST ${TITLE_LIST} --THRESHOLD_DEPTH ${THRESHOLD_DEPTH} --WES_TUMOR_BED ${WES_TUMOR_BED} --WES_DURA_VCF ${WES_DURA_VCF}


# WES_DURA VCF와  SINGLE_CELL_TUMOR_VCF 합쳐주기
echo -e "bcftools merge --force-samples  -l ${BCFTOOLS_MERGE_TXT} -Oz -o ${BCFTOOLS_MERGE_OUTPUT_VCF_GZ}"
bcftools merge --force-samples -l ${BCFTOOLS_MERGE_TXT} -Oz -o ${BCFTOOLS_MERGE_OUTPUT_VCF_GZ}
gunzip -c -f ${BCFTOOLS_MERGE_OUTPUT_VCF_GZ} >  ${BCFTOOLS_MERGE_OUTPUT_VCF_GZ%.gz}

# Heatmap 그리기
echo -e "python3 "06.HC_pipe_05.draw_integrated_heatmap.py" --BCFTOOLS_MERGE_OUTPUT_VCF  ${BCFTOOLS_MERGE_OUTPUT_VCF_GZ%.gz}  --OUTPUT_HEATMAP_PATH ${OUTPUT_HEATMAP_PATH}"
python3 "06.HC_pipe_05.draw_integrated_heatmap.py" --BCFTOOLS_MERGE_OUTPUT_VCF  ${BCFTOOLS_MERGE_OUTPUT_VCF_GZ%.gz}  --OUTPUT_HEATMAP_PATH ${OUTPUT_HEATMAP_PATH} 

# python3 "draw_integrated_heatmap.py" --BCFTOOLS_MERGE_OUTPUT_VCF /data/project/Meningioma/61.Lowinput/02.PTA/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/230405.merge.vcf --OUTPUT_HEATMAP_PATH /data/project/Meningioma/61.Lowinput/02.PTA/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/230405.heatmap.pdf