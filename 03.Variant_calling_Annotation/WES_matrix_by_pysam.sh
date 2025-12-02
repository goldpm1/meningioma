#!/bin/bash
#$ -cwd
#$ -S /bin/bash


source ~/.bashrc
conda activate cnvpytor



#######################  1.  WES-Tumor   #######################

python3 "/data/project/Meningioma/script/03.Variant_calling_Annotation/WES_matrix_by_pysam.py" \
    --TISSUE "Tumor" \
    --BAMTYPE "WES" \
    --BAM_DIR "/data/project/Meningioma/02.Align/hg38/Tumor/05.Final_bam" \
    --EXCLUDE "230526_WGS,230323_1,230405_1,230802_1,230802_3,240403_1,240403_2,240403_4,240412_1,240412_3,240412_4,231101,240110,230812,240612,240911" \
    --INTERVAL "/data/project/Meningioma/91.Lollipop/driver_information.txt" \
    --VAF_DF "/data/project/Meningioma/07.pysam/Tumor_vaf_df.tsv" \
    --COUNT_DF "/data/project/Meningioma/07.pysam/Tumor_count_df.tsv" \
    --HEATMAP_PATH "/data/project/Meningioma/07.pysam/Tumor_heatmap.pdf" \
    --BAMSNAP_DIR "/data/project/Meningioma/07.pysam/BAMSNAP_Tumor" \
    --THRESHOLD_END 4 \
    --THRESHOLD_BQ 15 \
    --THRESHOLD_MULTIALLELIC 1 \
    --PT_SORT "True" \
    --SHOW_BAMSNAP "False"
> "/data/project/Meningioma/07.pysam/Tumor_WES.log"


# #######################  2.  WES-Dura   #######################
python3 "/data/project/Meningioma/script/03.Variant_calling_Annotation/WES_matrix_by_pysam.py" \
    --TISSUE "Dura" \
    --BAMTYPE "WES" \
    --BAM_DIR "/data/project/Meningioma/02.Align/hg38/Dura/05.Final_bam" \
    --EXCLUDE "230323_1,230405_1,230802_1,230802_3,240403_1,240403_2,240403_4,240412_1,240412_3,240412_4,231101,240110" \
    --INTERVAL "/data/project/Meningioma/91.Lollipop/driver_information.txt" \
    --VAF_DF "/data/project/Meningioma/07.pysam/Dura_vaf_df.tsv" \
    --COUNT_DF "/data/project/Meningioma/07.pysam/Dura_count_df.tsv" \
    --HEATMAP_PATH "/data/project/Meningioma/07.pysam/Dura_heatmap.pdf" \
    --BAMSNAP_DIR "/data/project/Meningioma/07.pysam/BAMSNAP_Dura"     \
    --THRESHOLD_END 4 \
    --THRESHOLD_BQ 15 \
    --THRESHOLD_MULTIALLELIC 1 \
    --PT_SORT "True" \
    --SHOW_BAMSNAP "False"
> "/data/project/Meningioma/07.pysam/Dura_WES.log"