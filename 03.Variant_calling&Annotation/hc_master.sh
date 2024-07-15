#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/02.Align"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "31.hc_call" "32.hc_VQSR" "33.hc_HF" "34.hc_vep"; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

hg="hg38"
REF="/home/goldpm1/reference/genome.fa"
#REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz"
    
HC_PATH="/data/project/Meningioma/06.hc"


sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬

for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102, 230127, 230323, 230419

    for TISSUE in Blood Tumor Dura Ventricle Cortex; do   # 각각에 대해서 다 돌려보고 germline concordance도 확인해보자
        echo -e ${Sample_ID}"_"${TISSUE}
        BAM_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE}"/05.Final_bam/"${Sample_ID}"_"$TISSUE".bam"

        if [ -f ${BAM_PATH} ]; then     # File이 있어야만 진행
            #01. HC call
            OUTPUT_HC=${HC_PATH}"/01.call/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".vcf"
            if [ ! -d ${OUTPUT_HC%/*} ] ; then
                mkdir -p ${OUTPUT_HC%/*}
            fi
            qsub -pe smp 6 -e $logPath"/31.hc_call" -o $logPath"/31.hc_call" -N 'hc_31.'${Sample_ID}"_"${TISSUE}  -hold_jid "doc_"${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/"hc_pipe_01.call.sh \
                --BAM_PATH ${BAM_PATH} --INTERVAL ${INTERVAL} --REF ${REF} --dbSNP ${dbSNP} --OUTPUT_HC ${OUTPUT_HC}

            #02. VQSR
            INPUT_VCF=${OUTPUT_HC}
            RECAL_FILE=${HC_PATH}"/02.VQSR/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".recal"
            TRANCHES_FILE=${HC_PATH}"/02.VQSR/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".tranches"
            OUTPUT_VCF_GZ=${HC_PATH}"/02.VQSR/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".vcf.gz"
            OUTPUT_VCF=${HC_PATH}"/02.VQSR/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".vcf"
            if [ ! -d ${OUTPUT_VCF%/*} ] ; then
            mkdir -p ${OUTPUT_VCF%/*}
            fi
            # qsub -pe smp 5 -e $logPath"/32.hc_VQSR" -o $logPath"/32.hc_VQSR" -N 'hc_32.'${Sample_ID}"_"${TISSUE} -hold_jid 'hc_31.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/hc_pipe_02.gvcf.sh" \
            #     --INPUT_VCF ${INPUT_VCF} --RECAL_FILE ${RECAL_FILE} --TRANCHES_FILE ${TRANCHES_FILE} --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} --OUTPUT_VCF ${OUTPUT_VCF} --REF ${REF} 

            #03. HF (DP > 100)
            INPUT_VCF_GZ=${HC_PATH}"/01.call/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".vcf.gz"
            OUTPUT_VCF_GZ=${HC_PATH}"/03.HF/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".DP100.vcf.gz"
            OUTPUT_VCF=${HC_PATH}"/03.HF/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".DP100.vcf"
            if [ ! -d ${OUTPUT_VCF%/*} ] ; then
            mkdir -p ${OUTPUT_VCF%/*}
            fi
            qsub -pe smp 1 -e $logPath"/33.hc_HF" -o $logPath"/33.hc_HF" -N 'hc_33.'${Sample_ID}"_"${TISSUE} -hold_jid 'hc_31.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/hc_pipe_03.HF.sh" \
                --INPUT_VCF_GZ ${INPUT_VCF_GZ}  --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} --OUTPUT_VCF ${OUTPUT_VCF} --REF ${REF} 

            #04. VEP
            INPUT_VCF=${HC_PATH}"/03.HF/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".DP100.vcf"
            OUTPUT_VCF=${HC_PATH}"/04.vep/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".DP100.vep.vcf"
            if [ ! -d ${OUTPUT_VCF%/*} ] ; then
            mkdir -p ${OUTPUT_VCF%/*}
            fi
            qsub -pe smp 6 -e $logPath"/34.hc_vep" -o $logPath"/34.hc_vep" -N 'hc_34.'${Sample_ID}"_"${TISSUE} -hold_jid 'hc_33.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/hc_pipe_04.vep.sh" \
            --REF ${REF} --INPUT_VCF ${INPUT_VCF} --OUTPUT_VCF ${OUTPUT_VCF}
        fi
    done
done