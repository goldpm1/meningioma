#!/bin/bash
#$ -S /bin/bash
#$ -cwd

REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
REF_hg19="/home/goldpm1/reference/Broadhg19/Homo_sapiens_assembly19.fasta"             # 왜그런지는 모르겠지만 이 reference를 써야 돌아간다
REF_hg38="/home/goldpm1/reference/genome.fa"
hg="hg38"

DIR="/data/project/Meningioma"
PYCLONEVI_DIR="/data/project/Meningioma/31.Clonality"

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"




if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "pycl_12.pyclonevi_tumoronly" "pycl_13.vis_tumoronly" ; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done



source /home/goldpm1/.bashrc
conda deactivate

#################################################### Master Node (conda deactivate)  ##############################################

for Sample_ID in 190426 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT 241127 241127_FT 241127_FTRt 241127_Para \
    190524 221026 221102 221202 230127 230303_1 230323_2 230323_11 230419 230526 230802_1 230802_2 230802_3 230812 230920 231006 231025 231101 231206  240110 240325 240403_1 240403_2 240403_3 240412_1 240412_2 240412_3 240612 240911 241023 241023_1 \
    220930_1 220930_2 220930_3 220930_4  230405_2 230822 240320 241023_1 \
    250207 250207_Ant  \
    241016 241211 250425 250212 250502 250509 ; do 

    FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH=${PYCLONEVI_DIR}"/11.make_matrix_tumoronly/"${Sample_ID}"/"${Sample_ID}".facetcnv_to_pyclonevi.tsv"     # Sample_ID (Tumor) , TISSUE (Dura) 를 각각 의미
    FACETCNV_TO_PYCLONEVI_OUTPUT_PATH=${PYCLONEVI_DIR}"/12.pyclonevi_tumoronly/"${Sample_ID}"/"${Sample_ID}".facetcnv_to_pyclonevi.tsv"

    if [[ -f ${FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH} ]]; then

        ############ 02. pyclonevi 돌리기 ############
        if [ ! -d ${PYCLONEVI_DIR}"/12.pyclonevi_tumoronly/"${Sample_ID} ] ; then
            mkdir -p ${PYCLONEVI_DIR}"/12.pyclonevi_tumoronly/"${Sample_ID}
        fi

        qsub -pe smp 3  -o $logPath"/pycl_12.pyclonevi_tumoronly" -e $logPath"/pycl_12.pyclonevi_tumoronly" -N "pycl_12.fac_"${Sample_ID} ${CURRENT_PATH}"/pyclonevi_pipe_02.pyclonevi.sh"  \
            --INPUT_TSV ${FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH} \
            --OUTPUT_H5 ${PYCLONEVI_DIR}"/12.pyclonevi_tumoronly/"${Sample_ID}"/"${Sample_ID}".facetcnv_to_pyclonevi.h5" \
            --OUTPUT_TSV ${FACETCNV_TO_PYCLONEVI_OUTPUT_PATH}


        #############  13. visualization (CNV를 고려하기는 힘들다.)  ########################
        OUTPUT_VIS_FIG1=${PYCLONEVI_DIR}"/12.pyclonevi_tumoronly/"${Sample_ID}"/"${Sample_ID}".visualization.pdf"
        OUTPUT_VIS_FIG2=${PYCLONEVI_DIR}"/12.pyclonevi_tumoronly/"${Sample_ID}"/"${Sample_ID}".visualization-scaled.pdf"
        OUTPUT_VIS_DF=${PYCLONEVI_DIR}"/12.pyclonevi_tumoronly/"${Sample_ID}"/"${Sample_ID}".visualization_df.tsv"
        OUTPUT_DIR1=${PYCLONEVI_DIR}"/12.pyclonevi_tumoronly/decomposed"
        OUTPUT_DIR2=${PYCLONEVI_DIR}"/12.pyclonevi_tumoronly/scaled"

        for folder in ${OUTPUT_DIR1} ${OUTPUT_DIR2} ; do
            if [ ! -d ${folder} ] ; then
                mkdir -p ${folder}
            fi
        done

        qsub -pe smp 3  -o $logPath"/pycl_13.vis_tumoronly" -e $logPath"/pycl_13.vis_tumoronly" -N "pycl_13.vis_"${Sample_ID} -hold_jid "pycl_12.fac_"${Sample_ID} ${CURRENT_PATH}"/pyclonevi_pipe_13.visualization_tumoronly.sh"  \
            --Sample_ID ${Sample_ID} \
            --FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH} \
            --FACETCNV_TO_PYCLONEVI_OUTPUT_PATH ${FACETCNV_TO_PYCLONEVI_OUTPUT_PATH} \
            --OUTPUT_VIS_FIG1 ${OUTPUT_VIS_FIG1} \
            --OUTPUT_VIS_FIG2 ${OUTPUT_VIS_FIG2} \
            --OUTPUT_VIS_DF ${OUTPUT_VIS_DF} \
            --OUTPUT_DIR1 ${OUTPUT_DIR1}  \
            --OUTPUT_DIR2 ${OUTPUT_DIR2}

    fi


done
