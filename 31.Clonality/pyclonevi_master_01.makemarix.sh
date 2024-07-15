#!/bin/bash
#$ -S /bin/bash
#$ -cwd

REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
REF_hg19="/home/goldpm1/reference/Broadhg19/Homo_sapiens_assembly19.fasta"             # 왜그런지는 모르겠지만 이 reference를 써야 돌아간다
REF_hg38="/home/goldpm1/reference/genome.fa"
hg="hg38"

PROJECT_DIR="/data/project/Meningioma"
CASE_BAMpath=${PROJECT_DIR}"/02.bam/case"
CONTROL_BAMpath=${PROJECT_DIR}"/02.bam/control"

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

BAM_DIR="/data/project/Meningioma/02.Align"
MUTECT_DIR="/data/project/Meningioma/04.mutect"
HC_DIR="/data/project/Meningioma/06.hc"
GVCF_DIR="/data/project/Meningioma/05.gvcf/02.remove_nonref"
PYCLONEVI_DIR="/data/project/Meningioma/31.Clonality"
FACETCNV_DIR="/data/project/Meningioma/11.cnv/5.facetcnv"
SEQUENZA_DIR="/data/project/Meningioma/11.cnv/2.sequenza"


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


##################################################### Computer Node ##############################################


#for Sample_ID in 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT 220930 221026 221102 221202 230127 230303 230323_2 230405_2 230419 230526 230822 230920; do
for Sample_ID in 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT; do
    ## 01.Blood에서 30개의 HC call을 뽑아서 bed file꼴로 출력하기  (이젠 이건 필요없을듯)
    RANDOM_PICK="10"
    HC_BLOOD_RANDOM_PICK_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}".HC.random_pick_"${RANDOM_PICK}".bed"
    # python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.selectHCfromBlood.py" \
    #     --RANDOM_PICK ${RANDOM_PICK} \
    #     --HC_OUTPUT_PATH ${HC_DIR}"/03.HF/"${Sample_ID}"/Blood/"${Sample_ID}"_Blood.DP100.vcf" \
    #     --HC_BLOOD_RANDOM_PICK_PATH ${HC_BLOOD_RANDOM_PICK_PATH}



    # 02. TSV 파일 만들기 (Dura, Tumor에 관계없이 하나에 다 합치는 작업)
    SEQUENZA_TO_PYCLONEVI_MATRIX_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}".sequenza_to_pyclonevi.tsv"
    FACETCNV_TO_PYCLONEVI_MATRIX_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}".facetcnv_to_pyclonevi.tsv"
    if [ ! -d ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH%/*} ] ; then
        mkdir -p ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH%/*}
    fi
    rm -rf ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH}   # 지워야 계속 add 한다

    echo -e ${Sample_ID}
    for TISSUE in Tumor Dura; do   #Tumor Dura  여러 샘플 있어도 계속 쌓아도 됨
        SEQUENZA_MUTATION_PATH=${SEQUENZA_DIR}"/hg19to38/"${Sample_ID}"_"${TISSUE}"_mutations.txt"
        SEQUENZA_SEGMENT_PATH=${SEQUENZA_DIR}"/hg19to38/"${Sample_ID}"_"${TISSUE}"_segments.txt"
        SEQUENZA_PLOIDY_PATH=${SEQUENZA_DIR}"/hg19/"${Sample_ID}"_"${TISSUE}"_confints_CP.txt"
        SEQUENZA_PURITY_PLOIDY_PATH=${SEQUENZA_DIR}"/hg19/"${Sample_ID}"_"${TISSUE}"_purity_ploidy.txt"
        FACETCNV_OUTPUT_PATH=${FACETCNV_DIR}"/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}".vcf.gz"
        FACETCNV_PURITY_PLODY_PATH=${FACETCNV_DIR}"/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_purity_ploidy.txt"

        FACETCNV_TO_BED_DF_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".facetcnv_to_bed_df.tsv"
        
        if [ "${TISSUE}" == "Tumor" ]; then          
            if [[ "${Sample_ID}" =~ "190426" ]]; then  #190426_FFT_Tumor
                MUTECT_OUTPUT_PATH=${MUTECT_DIR}"/02.PASS/190426_Tumor_"${Sample_ID#190426_}".MT2.FMC.HF.RMBLACK.vep.vcf"
                SEQUENZA_MUTATION_PATH=${SEQUENZA_DIR}"/hg19to38/"${Sample_ID%_*}"_"${TISSUE}"_"${Sample_ID#*_}"_mutations.txt"
                SEQUENZA_SEGMENT_PATH=${SEQUENZA_DIR}"/hg19to38/"${Sample_ID%_*}"_"${TISSUE}"_"${Sample_ID#*_}"_segments.txt"
                SEQUENZA_PLOIDY_PATH=${SEQUENZA_DIR}"/hg19/"${Sample_ID%_*}"_"${TISSUE}"_"${Sample_ID#*_}"_confints_CP.txt"
                SEQUENZA_PURITY_PLOIDY_PATH=${SEQUENZA_DIR}"/hg19/"${Sample_ID%_*}"_"${TISSUE}"_"${Sample_ID#*_}"_purity_ploidy.txt"
                FACETCNV_OUTPUT_PATH=${FACETCNV_DIR}"/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}".vcf.gz"
                FACETCNV_PURITY_PLODY_PATH=${FACETCNV_DIR}"/"${Sample_ID%_*}"/"${TISSUE}"/"${Sample_ID}"_purity_ploidy.txt"

                FACETCNV_TO_BED_DF_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".facetcnv_to_bed_df.tsv"
            else
                MUTECT_OUTPUT_PATH=${MUTECT_DIR}"/02.PASS/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.vep.vcf"
            fi
        elif  [ "${TISSUE}" != "Tumor" ]; then
            if [[ "${Sample_ID}" =~ "190426" ]]; then        # 190426_Dura
                MUTECT_OUTPUT_PATH=${MUTECT_DIR}"/04.Other_rescue/190426_Dura.MT2.FMC.HF.RMBLACK.rescue.vep.vcf"
                SEQUENZA_MUTATION_PATH=${SEQUENZA_DIR}"/hg19to38/"${Sample_ID%_*}"_"${TISSUE}"_mutations.txt"
                SEQUENZA_SEGMENT_PATH=${SEQUENZA_DIR}"/hg19to38/"${Sample_ID%_*}"_"${TISSUE}"_segments.txt"
                SEQUENZA_PLOIDY_PATH=${SEQUENZA_DIR}"/hg19/"${Sample_ID%_*}"_"${TISSUE}"_confints_CP.txt"
                SEQUENZA_PURITY_PLOIDY_PATH=${SEQUENZA_DIR}"/hg19/"${Sample_ID%_*}"_"${TISSUE}"_purity_ploidy.txt"
                FACETCNV_OUTPUT_PATH=${FACETCNV_DIR}"/"${Sample_ID%_*}"/"${TISSUE}"/"${Sample_ID%_*}".vcf.gz"
                FACETCNV_PURITY_PLODY_PATH=${FACETCNV_DIR}"/"${Sample_ID%_*}"/"${TISSUE}"/"${Sample_ID%_*}"_purity_ploidy.txt"

                FACETCNV_TO_BED_DF_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".facetcnv_to_bed_df.tsv"
            else
                MUTECT_OUTPUT_PATH=${MUTECT_DIR}"/04.Other_rescue/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.vcf"
            fi
        fi
        HC_OUTPUT_PATH=${HC_DIR}"/04.vep/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"${TISSUE}".DP100.vep.vcf"



        python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.makematrix_sequenza_to_pyclonevi.py"  \
            --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} \
            --SEQUENZA_MUTATION_PATH ${SEQUENZA_MUTATION_PATH} --SEQUENZA_SEGMENT_PATH ${SEQUENZA_SEGMENT_PATH} --SEQUENZA_TO_PYCLONEVI_MATRIX_PATH ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} --SEQUENZA_PLOIDY_PATH ${SEQUENZA_PLOIDY_PATH}  --SEQUENZA_PURITY_PLOIDY_PATH ${SEQUENZA_PURITY_PLOIDY_PATH} \
            --MUTECT_OUTPUT_PATH ${MUTECT_OUTPUT_PATH} \
            --HC_OUTPUT_PATH ${HC_OUTPUT_PATH} --HC_BLOOD_RANDOM_PICK_PATH ${HC_BLOOD_RANDOM_PICK_PATH}
        # echo -e " python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.makematrix_sequenza_to_pyclonevi.py"  \
        #     --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} \
        #     --SEQUENZA_MUTATION_PATH ${SEQUENZA_MUTATION_PATH} --SEQUENZA_SEGMENT_PATH ${SEQUENZA_SEGMENT_PATH} --SEQUENZA_TO_PYCLONEVI_MATRIX_PATH ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} --SEQUENZA_PLOIDY_PATH ${SEQUENZA_PLOIDY_PATH}  --SEQUENZA_PURITY_PLOIDY_PATH ${SEQUENZA_PURITY_PLOIDY_PATH} \
        #     --MUTECT_OUTPUT_PATH ${MUTECT_OUTPUT_PATH} \
        #     --HC_OUTPUT_PATH ${HC_OUTPUT_PATH} --HC_BLOOD_RANDOM_PICK_PATH ${HC_BLOOD_RANDOM_PICK_PATH}"


        python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.makematrix_facetcnv_to_pyclonevi.py"  \
            --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} \
            --FACETCNV_OUTPUT_PATH ${FACETCNV_OUTPUT_PATH} --FACETCNV_TO_BED_DF_PATH ${FACETCNV_TO_BED_DF_PATH} --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} --FACETCNV_PURITY_PLODY_PATH ${FACETCNV_PURITY_PLODY_PATH} \
            --MUTECT_OUTPUT_PATH ${MUTECT_OUTPUT_PATH} \
            --HC_OUTPUT_PATH ${HC_OUTPUT_PATH} --HC_BLOOD_RANDOM_PICK_PATH ${HC_BLOOD_RANDOM_PICK_PATH}



    done

    # 03. 보기 좋게 Sort 하기 + PycloneVI를 위해 axis mutation (unique mutation) 도 살려줄지 결정하기
    python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.rescue+order.py" \
            --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} \
            --SEQUENZA_TO_PYCLONEVI_MATRIX_PATH ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} \
            --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
            --RESCUE_UNIQUEMUTATION True

done

