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


refine_log_dir() {
    local sublog="$1"
    local sublog_path="${logPath}/${sublog}"
    if [ -e "$sublog_path" ] ; then
        rm -rf "$sublog_path"
    fi
    if [ ! -d "$sublog_path" ] ; then
        mkdir -p "$sublog_path"
    fi
}


sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬


##################################################### Computer Node ##############################################


for Sample_ID in 190426 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT 241127 241127_FT 241127_FTRt 241127_Para  230323_2 230323_11 \
    190524  221102  230127  230419 230526 230802_1 230802_2 230802_3 230812 230822 230920 231025 231101 231206  240110 240320 240325 240403_1 240403_2 240403_3 240412_1 240412_2 240412_3 241023   \
    221202 231006 230303_1 230405_2 \
    220930_1 220930_2 220930_3 220930_4 221026 241016 241211 250425 250212 250502 250509 ; do 
    

    MUTECT_TUMOR_PATH=${MUTECT_DIR}"/02.PASS/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.vep.vcf"
    FACETCNV_TUMOR_OUTPUT_PATH=${FACETCNV_DIR}"/"${Sample_ID}"/Tumor/"${Sample_ID}".vcf.gz"
    FACETCNV_TUMOR_PURITY_PLODY_PATH=${FACETCNV_DIR}"/"${Sample_ID}"/Tumor/"${Sample_ID}"_purity_ploidy.txt"
    FACETCNV_TUMOR_TO_BED_DF_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}"_Tumor.facetcnv_to_bed_df.tsv"


    if [[ -f ${MUTECT_TUMOR_PATH} ]]; then

        for TISSUE in Dura DT Falx Falx_1 Falx_2 Dura_Ant; do     # Tumor x Dura 조합을 보기 위해
            MUTECT_DURA_PATH=${MUTECT_DIR}"/04.Other_rescue/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.vcf"
            FACETCNV_DURA_OUTPUT_PATH=${FACETCNV_DIR}"/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}".vcf.gz"
            FACETCNV_DURA_PURITY_PLODY_PATH=${FACETCNV_DIR}"/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_purity_ploidy.txt"
            FACETCNV_DURA_TO_BED_DF_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".facetcnv_to_bed_df.tsv"

            if [[ "${Sample_ID}" =~ "190426" ]]; then        # 190426_Dura
                #MUTECT_DURA_PATH=${MUTECT_DIR}"/04.Other_rescue/190426_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.vcf"
                FACETCNV_DURA_OUTPUT_PATH=${FACETCNV_DIR}"/"${Sample_ID%_*}"/"${TISSUE}"/"${Sample_ID%_*}".vcf.gz"
                FACETCNV_DURA_PURITY_PLODY_PATH=${FACETCNV_DIR}"/"${Sample_ID%_*}"/"${TISSUE}"/"${Sample_ID%_*}"_purity_ploidy.txt"
            elif [[ "${Sample_ID}" =~ "241127" ]]; then
                #MUTECT_DURA_PATH=${MUTECT_DIR}"/04.Other_rescue/241127_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.vcf"
                FACETCNV_DURA_OUTPUT_PATH=${FACETCNV_DIR}"/241127/"${TISSUE}"/241127.vcf.gz"
                FACETCNV_DURA_PURITY_PLODY_PATH=${FACETCNV_DIR}"/241127/"${TISSUE}"/241127_purity_ploidy.txt"
            elif [[ "${Sample_ID}" =~ "250207" ]]; then     # 250207
                #MUTECT_DURA_PATH=${MUTECT_DIR}"/04.Other_rescue/250207_Ant_"${TISSUE}".MT2.FMC.HF.RMBLACK.rescue.vep.vcf"
                FACETCNV_DURA_OUTPUT_PATH=${FACETCNV_DIR}"/250207/"${TISSUE}"/250207.vcf.gz"
                FACETCNV_DURA_PURITY_PLODY_PATH=${FACETCNV_DIR}"/250207/"${TISSUE}"/250207_purity_ploidy.txt"
            fi


            if [[ -f ${MUTECT_DURA_PATH} ]]; then
                echo -e "y axis : ${Sample_ID}_Tumor\t\tx axis : ${TISSUE}\t\t$(date)"

                # facetcnv_to_pyclonevi.tsv 파일 만들기
                FACETCNV_TO_PYCLONEVI_MATRIX_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".facetcnv_to_pyclonevi.tsv"     # Sample_ID (Tumor) , TISSUE (Dura) 를 각각 의미
                FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH=${PYCLONEVI_DIR}"/11.make_matrix_tumoronly/"${Sample_ID}"/"${Sample_ID}".facetcnv_to_pyclonevi.tsv"     # Sample_ID (Tumor) , TISSUE (Dura) 를 각각 의미
                if [ ! -d ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH%/*} ] ; then
                    mkdir -p ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH%/*}
                fi
                if [ ! -d ${FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH%/*} ] ; then
                    mkdir -p ${FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH%/*}
                fi
                # 빈 파일로 만들어줘야 계속 Tumor-dura 순으로 add 한다
                rm -rf ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH}  ${FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH}  
                touch ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} ${FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH}


                ###### Tumor 먼저 쌓아주기 (x축, open ("a") )  # ( ${Sample_ID}"_"${TISSUE}".facetcnv_to_pyclonevi.tsv" &   ${Sample_ID}"_Tumor.facetcnv_to_bed_df.tsv" )
                python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.makematrix_facetcnv_to_pyclonevi.py"  \
                    --Sample_ID ${Sample_ID} --TISSUE "Tumor" \
                    --FACETCNV_OUTPUT_PATH ${FACETCNV_TUMOR_OUTPUT_PATH} --FACETCNV_PURITY_PLODY_PATH ${FACETCNV_TUMOR_PURITY_PLODY_PATH}  --MUTECT_PATH ${MUTECT_TUMOR_PATH} \
                    --FACETCNV_TO_BED_DF_PATH ${FACETCNV_TUMOR_TO_BED_DF_PATH} --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
                    --REMOVE_TEMP True


                ###### Dura 쌓아주기 (y축)  # ( ${Sample_ID}"_"${TISSUE}".facetcnv_to_pyclonevi.tsv"  )
                python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.makematrix_facetcnv_to_pyclonevi.py"  \
                    --Sample_ID ${Sample_ID} --TISSUE "Dura" \
                    --FACETCNV_OUTPUT_PATH ${FACETCNV_DURA_OUTPUT_PATH} --FACETCNV_PURITY_PLODY_PATH ${FACETCNV_DURA_PURITY_PLODY_PATH}  --MUTECT_PATH ${MUTECT_DURA_PATH} \
                    --FACETCNV_TO_BED_DF_PATH ${FACETCNV_DURA_TO_BED_DF_PATH} --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
                    --REMOVE_TEMP True

                ###### 보기 좋게 Sort 하기 + PycloneVI를 위해 axis mutation (unique mutation) 도 살려줄지 결정하기  ( ${Sample_ID}"_"${TISSUE}".facetcnv_to_pyclonevi.tsv"  )
                python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.rescue+order.py" \
                    --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
                    --RESCUE_UNIQUEMUTATION "True"

                cp ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} ${FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH}
                python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.rescue+order.py" \
                    --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH} \
                    --RESCUE_UNIQUEMUTATION "False"
            fi

        done
    fi
done








#SEQUENZA_TO_PYCLONEVI_MATRIX_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}".sequenza_to_pyclonevi.tsv"
#HC_OUTPUT_PATH=${HC_DIR}"/04.vep/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"${TISSUE}".DP100.vep.vcf"

# ## 01.Blood에서 30개의 HC call을 뽑아서 bed file꼴로 출력하기  (이젠 이건 필요없을듯)
# RANDOM_PICK="10"
# HC_BLOOD_RANDOM_PICK_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}".HC.random_pick_"${RANDOM_PICK}".bed"
# # python3 ${CURRENT_PATH}"/pyclonevi_pipe_01.selectHCfromBlood.py" \
# #     --RANDOM_PICK ${RANDOM_PICK} \
# #     --HC_OUTPUT_PATH ${HC_DIR}"/03.HF/"${Sample_ID}"/Blood/"${Sample_ID}"_Blood.DP100.vcf" \
# #     --HC_BLOOD_RANDOM_PICK_PATH ${HC_BLOOD_RANDOM_PICK_PATH}


# # 예외상황
# if [[ "${Sample_ID}" =~ "190426" ]]; then  # 190426_FFT
#     FACETCNV_OUTPUT_PATH=${FACETCNV_DIR}"/"${Sample_ID%_*}"/Tumor_"${Sample_ID#190426_}"/"${Sample_ID%_*}".vcf.gz"
#     FACETCNV_PURITY_PLODY_PATH=${FACETCNV_DIR}"/"${Sample_ID%_*}"/Tumor_"${Sample_ID#190426_}"/"${Sample_ID%_*}"_purity_ploidy.txt"
#     FACETCNV_TO_BED_DF_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE_BROAD}".facetcnv_to_bed_df.tsv"
# fi
# if [[ "${Sample_ID}" =~ "241127_" ]]; then  #241127_FT
#     FACETCNV_OUTPUT_PATH=${FACETCNV_DIR}"/241127/Tumor/241127.vcf.gz"
#     FACETCNV_PURITY_PLODY_PATH=${FACETCNV_DIR}"/241127/Tumor/241127_purity_ploidy.txt"
#     FACETCNV_TO_BED_DF_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE_BROAD}".facetcnv_to_bed_df.tsv"
# fi