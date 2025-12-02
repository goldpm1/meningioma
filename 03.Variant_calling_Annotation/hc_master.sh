#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

BAM_DIR="/data/project/Meningioma/02.Align"
HC_DIR="/data/project/Meningioma/06.hc"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "31.hc_call" "32.hc_VQSR" "33.hc_HF" "34.hc_vep" "35.liftover"; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

hg="hg38"
if [ ${hg} == "hg38" ]; then
    REF="/home/goldpm1/reference/genome.fa"
    #REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
    INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
    INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_V8_hg38/S33613367_Covered.bed"      # SureSelect V8
    dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz"
elif [ ${hg} == "hg19" ]; then
    REF="/home/goldpm1/reference/hg19/hg19.fa"
    INTERVAL="/home/goldpm1/resources/SureSelectXTHumanAllExonV5Kit_hg19/S04380110_Covered.bed"
    dbSNP="/data/public/dbSNP/b154/GRCh37/GCF_000001405.25.re.common.vcf.gz"
fi


#for Sample_ID in 250310 250319 250407 250408 250428 250509 250513 250515 250520 250526_Left 250526_Right 250527 250530  250602 250605 250609  ; do 
for Sample_ID in 220930_3 221021 231011 231208 231220 240112 250227 250617 250627 250704 250716 250723_Ant 250723_Post 250725; do 
    #for TISSUE in Tumor Tumor_FFT Tumor_PCT Tumor_PFT Tumor_PP Tumor_PT Dura Falx Ventricle Cortex Tumor_FT Tumor_FTRt Tumor_Para; do   #Tumor Dura
    for TISSUE in Tumor ; do   #Tumor Dura
        if [[ $TISSUE == *"Falx"* ]]; then
            TISSUE_BROAD="Dura"
        else
            TISSUE_BROAD=${TISSUE%%_*}
        fi

        BAM_PATH=${BAM_DIR}"/"${hg}"/"${TISSUE_BROAD}"/05.Final_bam/"${Sample_ID}"_"$TISSUE".bam"

        if [ -f ${BAM_PATH} ]; then     # File이 있어야만 진행

            #echo -e ${Sample_ID}"_"${TISSUE}

            #01. HC call
            OUTPUT_HC=${HC_DIR}"/"${hg}"/01.call/"${Sample_ID}"/"${TISSUE_BROAD}"/"${Sample_ID}"_"$TISSUE".vcf"
            if [ ! -d ${OUTPUT_HC%/*} ] ; then
                mkdir -p ${OUTPUT_HC%/*}
            fi
            qsub -pe smp 6 -e $logPath"/31.hc_call" -o $logPath"/31.hc_call" -N 'hc_31.'${Sample_ID}"_"${TISSUE}  -hold_jid "doc_"${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/"hc_pipe_01.call.sh \
                --BAM_PATH ${BAM_PATH} --INTERVAL ${INTERVAL} --REF ${REF} --dbSNP ${dbSNP} --OUTPUT_HC ${OUTPUT_HC}

            #02. VQSR
            INPUT_VCF=${OUTPUT_HC}
            RECAL_FILE=${HC_DIR}"/"${hg}"/02.VQSR/"${Sample_ID}"/"${TISSUE_BROAD}"/"${Sample_ID}"_"$TISSUE".recal"
            TRANCHES_FILE=${HC_DIR}"/"${hg}"/02.VQSR/"${Sample_ID}"/"${TISSUE_BROAD}"/"${Sample_ID}"_"$TISSUE".tranches"
            OUTPUT_VCF_GZ=${HC_DIR}"/"${hg}"/02.VQSR/"${Sample_ID}"/"${TISSUE_BROAD}"/"${Sample_ID}"_"$TISSUE".vcf.gz"
            OUTPUT_VCF=${HC_DIR}"/"${hg}"/02.VQSR/"${Sample_ID}"/"${TISSUE_BROAD}"/"${Sample_ID}"_"$TISSUE".vcf"
            if [ ! -d ${OUTPUT_VCF%/*} ] ; then
                mkdir -p ${OUTPUT_VCF%/*}
            fi
            qsub -pe smp 5 -e $logPath"/32.hc_VQSR" -o $logPath"/32.hc_VQSR" -N 'hc_32.'${Sample_ID}"_"${TISSUE} -hold_jid 'hc_31.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/hc_pipe_02.gvcf.sh" \
                --INPUT_VCF ${INPUT_VCF} --RECAL_FILE ${RECAL_FILE} --TRANCHES_FILE ${TRANCHES_FILE} --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} --OUTPUT_VCF ${OUTPUT_VCF} --REF ${REF} 

            #03. HF (DP > 100)
            INPUT_VCF_GZ=${HC_DIR}"/"${hg}"/01.call/"${Sample_ID}"/"${TISSUE_BROAD}"/"${Sample_ID}"_"$TISSUE".vcf.gz"
            OUTPUT_VCF_GZ=${HC_DIR}"/"${hg}"/03.HF/"${Sample_ID}"/"${TISSUE_BROAD}"/"${Sample_ID}"_"$TISSUE".DP100.vcf.gz"
            OUTPUT_VCF=${HC_DIR}"/"${hg}"/03.HF/"${Sample_ID}"/"${TISSUE_BROAD}"/"${Sample_ID}"_"$TISSUE".DP100.vcf"
            if [ ! -d ${OUTPUT_VCF%/*} ] ; then
                mkdir -p ${OUTPUT_VCF%/*}
            fi
            qsub -pe smp 1 -e $logPath"/33.hc_HF" -o $logPath"/33.hc_HF" -N 'hc_33.'${Sample_ID}"_"${TISSUE} -hold_jid 'hc_31.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/hc_pipe_03.HF.sh" \
                --INPUT_VCF_GZ ${INPUT_VCF_GZ}  --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} --OUTPUT_VCF ${OUTPUT_VCF} --REF ${REF} 


            #04. VEP
            INPUT_VCF=${HC_DIR}"/"${hg}"/03.HF/"${Sample_ID}"/"${TISSUE_BROAD}"/"${Sample_ID}"_"$TISSUE_BROAD".DP100.vcf"
            OUTPUT_VCF=${HC_DIR}"/"${hg}"/04.vep/"${Sample_ID}"/"${TISSUE_BROAD}"/"${Sample_ID}"_"$TISSUE".DP100.vep.vcf"
            if [ ! -d ${OUTPUT_VCF%/*} ] ; then
                mkdir -p ${OUTPUT_VCF%/*}
            fi
            qsub -pe smp 6 -e $logPath"/34.hc_vep" -o $logPath"/34.hc_vep" -N 'hc_34.'${Sample_ID}"_"${TISSUE} -hold_jid 'hc_33.'${Sample_ID}"_"${TISSUE}  "/data/project/BrainTumor/91.script_YS/01.WES_basic/03.VariantCalling/vep_hg38.sh" \
                --INPUT_VCF ${INPUT_VCF} --OUTPUT_VCF ${OUTPUT_VCF}

            #05. hg38 -> hg19 liftover  (REF : Target REF)
            qsub -pe smp 2 -e ${logPath}"/35.liftover" -o ${logPath}"/35.liftover" -N  'hc_35.'${Sample_ID}"_"${TISSUE} -hold_jid  'hc_34.'${Sample_ID}"_"${TISSUE} "/data/project/BrainTumor/91.script_YS/99.Others/K1K_liftover_pipe.02.liftover.sh" \
                --INPUT_VCF ${HC_DIR}"/hg38/04.vep/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".DP100.vep.vcf"  \
                --OUTPUT_VCF ${HC_DIR}"/hg19/04.vep/"${Sample_ID}"/"${TISSUE}"/"${Sample_ID}"_"$TISSUE".DP100.vep.vcf"  \
                --REF "/home/goldpm1/reference/hg19/hg19.fa" \
                --LIFTOVER_CHAIN "/home/goldpm1/resources/hg38ToHg19.over.chain.gz"
        fi
    done
done