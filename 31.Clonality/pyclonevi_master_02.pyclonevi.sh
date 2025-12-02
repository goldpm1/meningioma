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





# source /home/goldpm1/.bashrc
# conda deactivate

#################################################### Master Node (conda deactivate)  ##############################################

refine_log_dir "pycl_02.pyclonevi" "pycl_03.vis"

for Sample_ID in 190426 190426_FFT 190426_PCT 190426_PFT 190426_PP 190426_PT 241127 241127_FT 241127_FTRt 241127_Para  230323_2 230323_11 \
    190524  221102  230127  230419 230526 230802_1 230802_2 230802_3 230812 230822 230920 231025 231101 231206  240110 240320 240325 240403_1 240403_2 240403_3 240412_1 240412_2 240412_3 241023   \
    221202 231006 230303_1 230405_2 \
    220930_1 220930_2 220930_3 220930_4 221026 241016 241211 250425 250212 250502 250509 ; do 
    
    for TISSUE in Dura DT Falx Falx_1 Falx_2 Dura_Ant; do     # Tumor x Dura 조합을 보기 위해

        FACETCNV_TO_PYCLONEVI_MATRIX_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".facetcnv_to_pyclonevi.tsv"
        FACETCNV_TO_PYCLONEVI_OUTPUT_PATH=${PYCLONEVI_DIR}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".facetcnv_to_pyclonevi.tsv"

        if [[ -f ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} ]]; then

            ############ 02. pyclonevi 돌리기 ############
            if [ ! -d ${PYCLONEVI_DIR}"/02.pyclonevi/"${Sample_ID} ] ; then
                mkdir -p ${PYCLONEVI_DIR}"/02.pyclonevi/"${Sample_ID}
            fi

            echo -e "y axis : ${Sample_ID}_Tumor\t\tx axis : ${TISSUE}\t\t$(date)"

            PYCL02_JOB=$(sbatch --job-name='pycl_02.fac_'${Sample_ID}'_'${TISSUE} --output=${logPath}"/pycl_02.pyclonevi/"${Sample_ID}"_"${TISSUE}".log" --error=${logPath}"/pycl_02.pyclonevi/"${Sample_ID}"_"${TISSUE}".err" \
                --cpus-per-task=3 --mem=14GB --time=0:30:00 --partition=cpu --qos nstumor \
                --wrap="bash ${CURRENT_PATH}/pyclonevi_pipe_02.pyclonevi.sh \
                    --INPUT_TSV ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
                    --OUTPUT_H5 ${PYCLONEVI_DIR}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".facetcnv_to_pyclonevi.h5" \
                    --OUTPUT_TSV ${FACETCNV_TO_PYCLONEVI_OUTPUT_PATH} ")
            PYCL02_JOBID=$(echo ${PYCL02_JOB} | awk '{print $NF}')


            #############  03. visualization (CNV를 고려하기는 힘들다.)  ########################
            OUTPUT_VIS_FIG1=${PYCLONEVI_DIR}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".visualization.pdf"
            OUTPUT_VIS_FIG2=${PYCLONEVI_DIR}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".visualization-scaled.pdf"
            OUTPUT_VIS_DF=${PYCLONEVI_DIR}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".visualization_df.tsv"
            OUTPUT_DIR1=${PYCLONEVI_DIR}"/02.pyclonevi/decomposed"
            OUTPUT_DIR2=${PYCLONEVI_DIR}"/02.pyclonevi/scaled"

            for folder in ${OUTPUT_DIR1} ${OUTPUT_DIR2} ; do
                if [ ! -d ${folder} ] ; then
                    mkdir -p ${folder}
                fi
            done

            PYCL03_JOB=$(sbatch  --dependency=afterok:${PYCL02_JOBID} --job-name='pycl_03.vis_'${Sample_ID}'_'${TISSUE} --output=${logPath}"/pycl_03.vis/"${Sample_ID}"_"${TISSUE}".log" --error=${logPath}"/pycl_03.vis/"${Sample_ID}"_"${TISSUE}".err" \
                --cpus-per-task=1 --mem=14GB --time=0:30:00 --partition=cpu --qos nstumor \
                --wrap="bash ${CURRENT_PATH}/pyclonevi_pipe_03.visualization.sh \
                    --Sample_ID ${Sample_ID}"_"${TISSUE} \
                    --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
                    --FACETCNV_TO_PYCLONEVI_OUTPUT_PATH ${FACETCNV_TO_PYCLONEVI_OUTPUT_PATH} \
                    --OUTPUT_VIS_FIG1 ${OUTPUT_VIS_FIG1} \
                    --OUTPUT_VIS_FIG2 ${OUTPUT_VIS_FIG2} \
                    --OUTPUT_VIS_DF ${OUTPUT_VIS_DF} \
                    --OUTPUT_DIR1 ${OUTPUT_DIR1}  \
                    --OUTPUT_DIR2 ${OUTPUT_DIR2} ")
            PYCL03_JOBID=$(echo ${PYCL03_JOB} | awk '{print $NF}')

        fi

    done

done
    
    
# SEQUENZA_TO_PYCLONEVI_MATRIX_PATH=${PYCLONEVI_DIR}"/01.make_matrix/"${Sample_ID}"/"${Sample_ID}".sequenza_to_pyclonevi.tsv"
# SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH=${PYCLONEVI_DIR}"/02.pyclonevi/"${Sample_ID}"/"${Sample_ID}".sequenza_to_pyclonevi.tsv"