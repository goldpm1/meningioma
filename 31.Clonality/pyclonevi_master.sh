#!/bin/bash
#$ -S /bin/bash
#$ -cwd

REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
REF_hg19="/home/goldpm1/reference/Broadhg19/Homo_sapiens_assembly19.fasta"             # 왜그런지는 모르겠지만 이 reference를 써야 돌아간다
REF_hg38="/home/goldpm1/reference/genome.fa"
hg="hg38"

DIR="/home/goldpm1/Meningioma"
CASE_BAMpath=$DIR"/02.bam/case"
CONTROL_BAMpath=$DIR"/02.bam/control"



CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/home/goldpm1/Meningioma/02.Align"
MUTECT_PATH="/home/goldpm1/Meningioma/04.mutect/04.rescue"
HC_PATH="/home/goldpm1/Meningioma/06.hc"
GVCF_PATH="/home/goldpm1/Meningioma/05.gvcf/02.remove_nonref"
PYCLONEVIpath="/home/goldpm1/Meningioma/31.Clonality"
FACETCNV_PATH="/home/goldpm1/Meningioma/11.cnv/5.facetcnv"
SEQUENZA_PATH="/home/goldpm1/Meningioma/11.cnv"


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


sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬


for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102


    # 01.Blood에서 30개의 HC call을 뽑아서 bed file꼴로 출력하기
    RANDOM_PICK="10"
    HC_BLOOD_RANDOM_PICK_PATH=${PYCLONEVIpath}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}".HC.random_pick_"${RANDOM_PICK}".bed"
    python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.selectHCfromBlood.py" \
        --RANDOM_PICK ${RANDOM_PICK} \
        --HC_OUTPUT_PATH ${HC_PATH}"/03.HF/"${Sample_ID}"/Blood/"${Sample_ID}"_Blood.DP100.vcf" \
        --HC_BLOOD_RANDOM_PICK_PATH ${HC_BLOOD_RANDOM_PICK_PATH}



    # 02. TSV 파일 만들기 (Dura, Tumor에 관계없이 하나에 다 합치는 작업)
    SEQUENZA_TO_PYCLONEVI_MATRIX_PATH=${PYCLONEVIpath}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}".sequenza_to_pyclonevi.tsv"
    FACETCNV_TO_PYCLONEVI_MATRIX_PATH=${PYCLONEVIpath}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}".facetcnv_to_pyclonevi.tsv"
    if [ ! -d ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH%/*} ] ; then
        mkdir -p ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH%/*}
    fi
    rm -rf ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH}   # 지워야 계속 add 한다

    for TISSUE in Tumor Dura; do   #Tumor Dura  여러 샘플 있어도 계속 쌓아도 됨
        echo -e ${Sample_ID}"_"${TISSUE}
        SEQUENZA_MUTATION_PATH=${SEQUENZA_PATH}"/2.sequenza/hg19to38/"${Sample_ID}"_"${TISSUE}"_mutations.txt"
        SEQUENZA_SEGMENT_PATH=${SEQUENZA_PATH}"/2.sequenza/hg19to38/"${Sample_ID}"_"${TISSUE}"_segments.txt"
        SEQUENZA_PLOIDY_PATH=${SEQUENZA_PATH}"/2.sequenza/hg19/"${Sample_ID}"_"${TISSUE}"_confints_CP.txt"
        SEQUENZA_PURITY_PLOIDY_PATH=${SEQUENZA_PATH}"/2.sequenza/hg19/"${Sample_ID}"_"${TISSUE}"_purity_ploidy.txt"
        FACETCNV_OUTPUT_PATH=${FACETCNV_PATH}"/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}".vcf.gz"
        FACETCNV_PURITY_PLODY_PATH=${FACETCNV_PATH}"/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_purity_ploidy.txt"
        FACETCNV_TO_BED_DF_PATH=${PYCLONEVIpath}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".facetcnv_to_bed_df.tsv"
        MUTECT_OUTPUT_PATH=${MUTECT_PATH}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.rescue.vcf"
        HC_OUTPUT_PATH=${HC_PATH}"/04.vep/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"${TISSUE}".DP100.vep.vcf"


        python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.makematrix_sequenza_to_pyclonevi.py"  \
            --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} \
            --SEQUENZA_MUTATION_PATH ${SEQUENZA_MUTATION_PATH} --SEQUENZA_SEGMENT_PATH ${SEQUENZA_SEGMENT_PATH} --SEQUENZA_TO_PYCLONEVI_MATRIX_PATH ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} --SEQUENZA_PLOIDY_PATH ${SEQUENZA_PLOIDY_PATH}  --SEQUENZA_PURITY_PLOIDY_PATH ${SEQUENZA_PURITY_PLOIDY_PATH} \
            --MUTECT_OUTPUT_PATH ${MUTECT_OUTPUT_PATH} \
            --HC_OUTPUT_PATH ${HC_OUTPUT_PATH} --HC_BLOOD_RANDOM_PICK_PATH ${HC_BLOOD_RANDOM_PICK_PATH}

        python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.makematrix_facetcnv_to_pyclonevi.py"  \
            --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} \
            --FACETCNV_OUTPUT_PATH ${FACETCNV_OUTPUT_PATH} --FACETCNV_TO_BED_DF_PATH ${FACETCNV_TO_BED_DF_PATH} --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} --FACETCNV_PURITY_PLODY_PATH ${FACETCNV_PURITY_PLODY_PATH} \
            --MUTECT_OUTPUT_PATH ${MUTECT_OUTPUT_PATH} \
            --HC_OUTPUT_PATH ${HC_OUTPUT_PATH} --HC_BLOOD_RANDOM_PICK_PATH ${HC_BLOOD_RANDOM_PICK_PATH}

    done

    # 03. 보기 좋게 Sort 하기
    python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.order_by_type.py" \
            --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} \
            --SEQUENZA_TO_PYCLONEVI_MATRIX_PATH ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} \
            --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH}

    
    # 04. pyclonevi 돌리기
    if [ ! -d ${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID} ] ; then
        mkdir -p ${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}
    fi
    qsub -pe smp 3  -o $logPath"/pycl_02.pyclonevi" -e $logPath"/pycl_02.pyclonevi" -N "pycl2_seq_"${Sample_ID} ${CURRENT_PATH}"/pyclonevi_pipe_02.pyclonevi.sh"  \
        --INPUT_TSV ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} \
        --OUTPUT_H5 ${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".sequenza_to_pyclonevi.h5" \
        --OUTPUT_TSV ${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".sequenza_to_pyclonevi.tsv"

    qsub -pe smp 3  -o $logPath"/pycl_02.pyclonevi" -e $logPath"/pycl_02.pyclonevi" -N "pycl2_fac_"${Sample_ID} ${CURRENT_PATH}"/pyclonevi_pipe_02.pyclonevi.sh"  \
        --INPUT_TSV ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
        --OUTPUT_H5 ${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".facetcnv_to_pyclonevi.h5" \
        --OUTPUT_TSV ${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".facetcnv_to_pyclonevi.tsv"


    #05. visualization (CNV를 고려하기는 힘들다.)
    OUTPUT_PATH_SHARED=${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".decomposed.png"
    OUTPUT_PATH_TOTAL=${PYCLONEVIpath}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".total.png"
    qsub -pe smp 3  -o $logPath"/pycl_03.vis" -e $logPath"/pycl_03.vis" -N "pycl3_vis_"${Sample_ID} -hold_jid "pycl2_seq_"${Sample_ID}",pycl2_fac_"${Sample_ID} ${CURRENT_PATH}"/pyclonevi_pipe_03.visualization.sh"  \
        --Sample_ID ${Sample_ID} \
        --SEQUENZA_TO_PYCLONEVI_MATRIX_PATH ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} \
        --SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH ${SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH} \
        --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
        --FACETCNV_TO_PYCLONEVI_OUTPUT_PATH ${FACETCNV_TO_PYCLONEVI_OUTPUT_PATH} \
        --OUTPUT_PATH_SHARED ${OUTPUT_PATH_SHARED} \
        --OUTPUT_PATH_TOTAL ${OUTPUT_PATH_TOTAL}

done
