#!/bin/bash
#$ -cwd
#$ -S /bin/bash


############### Germline Joint Diploid ##################


REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
REF_hg19="/home/goldpm1/reference/Broadhg19/Homo_sapiens_assembly19.fasta"             # 왜그런지는 모르겠지만 이 reference를 써야 돌아간다
#REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
REF_hg38="/home/goldpm1/reference/genome.fa"

INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed.gz"
INTERVAL="/home/goldpm1/resources/whole.exome.proteincoding.bed.gz"
dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz"

hg="hg38"

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

for sublog in "02.manta_multiple" ; do
    if [ -d $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

DATA_PATH="/data/project/Meningioma/02.Align"
MANTA_PATH="/data/project/Meningioma/21.SV/02.manta_multiple"

for subpath in "01.raw" "02.PASS"; do
    if [ ! -d $MANTA_PATH"/"$subpath ] ; then
        mkdir -p $MANTA_PATH"/"$subpath
    fi
done

sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬


for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102

    # Basic  Input
    CASE_BAM_PATH_TUMOR=${DATA_PATH}"/"${hg}"/Tumor/05.Final_bam/"${Sample_ID}"_Tumor.bam"
    CASE_BAM_PATH_DURA=${DATA_PATH}"/"${hg}"/Dura/05.Final_bam/"${Sample_ID}"_Dura.bam"
    CASE_BAM_PATH_VENTRICLE=${DATA_PATH}"/"${hg}"/Ventricle/05.Final_bam/"${Sample_ID}"_Ventricle.bam"
    if [ ! -f ${CASE_BAM_PATH_VENTRICLE} ]; then     # File이 있어야만 진행
        CASE_BAM_PATH_VENTRICLE="None"
    fi
    CASE_BAM_PATH_CORTEX=${DATA_PATH}"/"${hg}"/Cortex/05.Final_bam/"${Sample_ID}"_Cortex.bam"
    if [ ! -f ${CASE_BAM_PATH_CORTEX} ]; then     # File이 있어야만 진행
        CASE_BAM_PATH_CORTEX="None"
    fi
    
    CONTROL_BAM_PATH=${DATA_PATH}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"


    # Output
    OUTPUT_DIR=${MANTA_PATH}"/01.raw/"${Sample_ID}
    OUTPUT_PATH_SOMATIC=${MANTA_PATH}"/02.PASS/01.Somatic/"${Sample_ID}".Manta.Somatic.PASS.vcf"
    OUTPUT_PATH_SOMATIC_GZ=${OUTPUT_PATH_SOMATIC}".gz"
    OUTPUT_PATH_DIPLOID=${MANTA_PATH}"/02.PASS/02.Diploid/"${Sample_ID}".Manta.Diploid.PASS.vcf"
    OUTPUT_PATH_DIPLOID_GZ=${OUTPUT_PATH_DIPLOID}".gz"
    if [ -d ${OUTPUT_DIR} ] ; then
        rm -rf ${OUTPUT_DIR}
    fi
    for subpath in ${OUTPUT_DIR} ${OUTPUT_PATH_SOMATIC%/*} ${OUTPUT_PATH_DIPLOID%/*} ; do
        if [ ! -d $subpath ] ; then
            mkdir -p $subpath
        fi
    done

    ####################################### 3개만 있는 경우 (Tumor, Dura, Ventricle) #######################################
    if [[ ${CASE_BAM_PATH_VENTRICLE} != "None" && ${CASE_BAM_PATH_CORTEX} == "None" ]]; then
        echo -e "\n"${Sample_ID}" : Tumor, Dura, Ventricle"

        for NUM in "12" "13" "23"; do        # 12: Tumor-Dura,      13: Tumor-Ventricle         23: Dura-Ventricle
            if [[ ${NUM} == "12" ]]; then
                CASE_BAM_PATH1=${CASE_BAM_PATH_TUMOR}
                CASE_BAM_PATH2=${CASE_BAM_PATH_DURA}
            fi
            if [[ ${NUM} == "13" ]]; then
                CASE_BAM_PATH1=${CASE_BAM_PATH_TUMOR}
                CASE_BAM_PATH2=${CASE_BAM_PATH_VENTRICLE}
            fi
            if [[ ${NUM} == "23" ]]; then
                CASE_BAM_PATH1=${CASE_BAM_PATH_DURA}
                CASE_BAM_PATH2=${CASE_BAM_PATH_VENTRICLE}
            fi

            qsub -pe smp 3 -o $logPath"/02.manta_multiple" -e $logPath"/02.manta_multiple" -N "Man01_"${Sample_ID}  "manta_multiple_pipe_01.call.sh" \
                --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} --CASE_BAM_PATH1 ${CASE_BAM_PATH1}  --CASE_BAM_PATH2 ${CASE_BAM_PATH2} \
                --REF ${REF_hg38} --CALLREGIONS ${INTERVAL} --OUTPUT_DIR ${OUTPUT_DIR} \
                --OUTPUT_PATH_SOMATIC ${OUTPUT_PATH_SOMATIC} --OUTPUT_PATH_SOMATIC_GZ ${OUTPUT_PATH_SOMATIC_GZ} \
                --OUTPUT_PATH_DIPLOID ${OUTPUT_PATH_DIPLOID} --OUTPUT_PATH_DIPLOID_GZ ${OUTPUT_PATH_DIPLOID_GZ}

        done
    fi
    ####################################### 3개만 있는 경우 (Tumor, Dura, Cortex) #######################################
    if [[ ${CASE_BAM_PATH_CORTEX} != "None" && ${CASE_BAM_PATH_VENTRICLE} == "None" ]]; then
        echo -e "\n"${Sample_ID}" : Tumor, Dura, Cortex"

        for NUM in "12" "14" "24"; do        # 12: Tumor-Dura,      14: Tumor-Cortex         24: Dura-Cortex
            if [[ ${NUM} == "12" ]]; then
                CASE_BAM_PATH1=${CASE_BAM_PATH_TUMOR}
                CASE_BAM_PATH2=${CASE_BAM_PATH_DURA}
            fi
            if [[ ${NUM} == "14" ]]; then
                CASE_BAM_PATH1=${CASE_BAM_PATH_TUMOR}
                CASE_BAM_PATH2=${CASE_BAM_PATH_CORTEX}
            fi
            if [[ ${NUM} == "24" ]]; then
                CASE_BAM_PATH1=${CASE_BAM_PATH_DURA}
                CASE_BAM_PATH2=${CASE_BAM_PATH_CORTEX}
            fi

            qsub -pe smp 3 -o $logPath"/02.manta_multiple" -e $logPath"/02.manta_multiple" -N "Man01_"${Sample_ID}  "manta_multiple_pipe_01.call.sh" \
                --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} --CASE_BAM_PATH1 ${CASE_BAM_PATH1}  --CASE_BAM_PATH2 ${CASE_BAM_PATH2} \
                --REF ${REF_hg38} --CALLREGIONS ${INTERVAL} --OUTPUT_DIR ${OUTPUT_DIR} \
                --OUTPUT_PATH_SOMATIC ${OUTPUT_PATH_SOMATIC} --OUTPUT_PATH_SOMATIC_GZ ${OUTPUT_PATH_SOMATIC_GZ} \
                --OUTPUT_PATH_DIPLOID ${OUTPUT_PATH_DIPLOID} --OUTPUT_PATH_DIPLOID_GZ ${OUTPUT_PATH_DIPLOID_GZ}

        done
    fi

    # ####################################### 2개만 있는 경우 (Tumor, Dura) #######################################
    if [[ ${CASE_BAM_PATH_VENTRICLE} == "None" && ${CASE_BAM_PATH_CORTEX} == "None" ]]; then
        echo -e "\n"${Sample_ID}" : Tumor, Dura"

        for NUM in "12" ; do        # 12: Tumor-Dura
            if [[ ${NUM} == "12" ]]; then
                CASE_BAM_PATH1=${CASE_BAM_PATH_TUMOR}
                CASE_BAM_PATH2=${CASE_BAM_PATH_DURA}
            fi

            qsub -pe smp 3 -o $logPath"/02.manta_multiple" -e $logPath"/02.manta_multiple" -N "Man01_"${Sample_ID}  "manta_multiple_pipe_01.call.sh" \
                --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} --CASE_BAM_PATH1 ${CASE_BAM_PATH1}  --CASE_BAM_PATH2 ${CASE_BAM_PATH2} \
                --REF ${REF_hg38} --CALLREGIONS ${INTERVAL} --OUTPUT_DIR ${OUTPUT_DIR} \
                --OUTPUT_PATH_SOMATIC ${OUTPUT_PATH_SOMATIC} --OUTPUT_PATH_SOMATIC_GZ ${OUTPUT_PATH_SOMATIC_GZ} \
                --OUTPUT_PATH_DIPLOID ${OUTPUT_PATH_DIPLOID} --OUTPUT_PATH_DIPLOID_GZ ${OUTPUT_PATH_DIPLOID_GZ}
        done
    fi
    
done




# # 02. Pandas Dataframe으로 만들기

# for idx in ${!sample_name_LIST[@]}; do
#     Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102

#     INPUT_PATH_DIPLOID=${MANTA_PATH}"/02.PASS/02.Diploid/"${Sample_ID}".Manta.Diploid.PASS.vcf"
#     OUTPUT_PATH_DIPLOID=${MANTA_PATH}"/02.PASS/02.Diploid/"${Sample_ID}".Manta.Diploid.PASS.chr.vcf"
#     OUTPUT_INV_PATH_DIPLOID=${MANTA_PATH}"/02.PASS/02.Diploid/"${Sample_ID}".Manta.PASS.INV.vcf"
#     PANDAS_DIR=${MANTA_PATH}"/03.pandas/"${Sample_ID}

#     if [ -d ${PANDAS_DIR} ] ; then
#         rm -rf ${PANDAS_DIR}
#     fi
#     for subpath in ${PANDAS_DIR} ; do
#         if [ ! -d $subpath ] ; then
#             mkdir -p $subpath
#         fi
#     done


#     qsub -pe smp 1 -o $logPath"/02.manta_multiple" -e $logPath"/02.manta_multiple" -N "man02_"${Sample_ID}  "manta_multiple_pipe_02.somaticSV.sh" \
#         --INPUT_PATH ${INPUT_PATH_DIPLOID} --OUTPUT_PATH ${OUTPUT_PATH_DIPLOID} --OUTPUT_INV_PATH ${OUTPUT_INV_PATH_DIPLOID} --PANDAS_DIR ${PANDAS_DIR} \
#         --ID ${Sample_ID} --CONTROL_SAMPLE ${Sample_ID}"_Blood" --REF ${REF_hg38} 
# done