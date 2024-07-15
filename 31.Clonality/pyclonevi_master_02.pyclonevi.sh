#!/bin/bash
#$ -S /bin/bash
#$ -cwd

REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
REF_hg19="/home/goldpm1/reference/Broadhg19/Homo_sapiens_assembly19.fasta"             # 왜그런지는 모르겠지만 이 reference를 써야 돌아간다
REF_hg38="/home/goldpm1/reference/genome.fa"
hg="hg38"

DIR="/data/project/Meningioma"
CASE_BAMpath=$DIR"/02.bam/case"
CONTROL_BAMpath=$DIR"/02.bam/control"



CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/02.Align"
MUTECT_PATH="/data/project/Meningioma/04.mutect/04.rescue"
HC_PATH="/data/project/Meningioma/06.hc"
GVCF_PATH="/data/project/Meningioma/05.gvcf/02.remove_nonref"
PYCLONEVIpath="/data/project/Meningioma/31.Clonality"
FACETCNV_PATH="/data/project/Meningioma/11.cnv/5.facetcnv"
SEQUENZA_PATH="/data/project/Meningioma/11.cnv"


if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "pycl_01.make_matrix" "pycl_02.pyclonevi" "pycl_03.vis" ; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


# sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
# sample_name_LIST=(${sample_name_list// / })     # array로 만듬
# for idx in ${!sample_name_LIST[@]}; do
#     Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102



#################################################### Master Node ##############################################

for Sample_ID in 190426 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT 220930 221026 221102 221202 230127 230303 230323_2 230405_2 230419 230526 230822 230920; do
#for Sample_ID in 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT; do
    SEQUENZA_TO_PYCLONEVI_MATRIX_PATH=${PYCLONEVIpath}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}".sequenza_to_pyclonevi.tsv"
    FACETCNV_TO_PYCLONEVI_MATRIX_PATH=${PYCLONEVIpath}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}".facetcnv_to_pyclonevi.tsv"
    SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH=${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".sequenza_to_pyclonevi.tsv"
    FACETCNV_TO_PYCLONEVI_OUTPUT_PATH=${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".facetcnv_to_pyclonevi.tsv"


    # 04. pyclonevi 돌리기
    if [ ! -d ${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID} ] ; then
        mkdir -p ${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}
    fi

    # qsub -pe smp 3  -o $logPath"/pycl_02.pyclonevi" -e $logPath"/pycl_02.pyclonevi" -N "pycl_02.seq_"${Sample_ID} ${CURRENT_PATH}"/pyclonevi_pipe_02.pyclonevi.sh"  \
    #     --INPUT_TSV ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} \
    #     --OUTPUT_H5 ${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".sequenza_to_pyclonevi.h5" \
    #     --OUTPUT_TSV ${SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH}

    # qsub -pe smp 3  -o $logPath"/pycl_02.pyclonevi" -e $logPath"/pycl_02.pyclonevi" -N "pycl_02.fac_"${Sample_ID} ${CURRENT_PATH}"/pyclonevi_pipe_02.pyclonevi.sh"  \
    #     --INPUT_TSV ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
    #     --OUTPUT_H5 ${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".facetcnv_to_pyclonevi.h5" \
    #     --OUTPUT_TSV ${FACETCNV_TO_PYCLONEVI_OUTPUT_PATH}


    #05. visualization (CNV를 고려하기는 힘들다.)
    OUTPUT_PATH_SHARED=${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".decomposed.pdf"
    OUTPUT_DIR1=${PYCLONEVIpath}"/02.pyclonevi/decomposed"
    OUTPUT_DIR2=${PYCLONEVIpath}"/02.pyclonevi/scaled"

    if [ ! -d ${OUTPUT_DIR1} ] ; then
        mkdir ${OUTPUT_DIR1}
    fi
    if [ ! -d ${OUTPUT_DIR2} ] ; then
        mkdir ${OUTPUT_DIR2}
    fi

    qsub -pe smp 3  -o $logPath"/pycl_03.vis" -e $logPath"/pycl_03.vis" -N "pycl_03.vis_"${Sample_ID} -hold_jid "pycl_02.seq_"${Sample_ID}",pycl2_fac_"${Sample_ID} ${CURRENT_PATH}"/pyclonevi_pipe_03.visualization.sh"  \
        --Sample_ID ${Sample_ID} \
        --SEQUENZA_TO_PYCLONEVI_MATRIX_PATH ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} \
        --SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH ${SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH} \
        --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
        --FACETCNV_TO_PYCLONEVI_OUTPUT_PATH ${FACETCNV_TO_PYCLONEVI_OUTPUT_PATH} \
        --OUTPUT_PATH_SHARED ${OUTPUT_PATH_SHARED} \
        --OUTPUT_DIR1 ${OUTPUT_DIR1}  \
        --OUTPUT_DIR2 ${OUTPUT_DIR2}

done