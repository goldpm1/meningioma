#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/home/goldpm1/Meningioma/02.Align"
PON="/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
#REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
REF="/home/goldpm1/reference/genome.fa"
hg="hg38"
gnomad="/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
TMP_PATH=${DATA_PATH}"/temp"
INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi

for sublog in 01.MTcall 02.FMC_HF_RMBLACK 03.vep 04.rescue 06.maf 11.multipleMT 12.vep; do
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
    
    for TISSUE in Tumor Dura ; do   #Tumor Dura
        CASE_BAM_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
        CONTROL_BAM_PATH=${DATA_PATH}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
        OUTPUT_VCF_GZ=${DATA_PATH%/*}"/04.mutect/01.raw/"${Sample_ID}"_"${TISSUE}".vcf.gz"
        OUTPUT_FMC_PATH=${DATA_PATH%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.vcf"
        OUTPUT_FMC_HF_PATH=${DATA_PATH%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.vcf"
        OUTPUT_FMC_HF_RMBLACK_PATH=${DATA_PATH%/*}"/04.mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vcf"


        SAMPLE_THRESHOLD="Dura,Tumor"   # "all"
        DP_THRESHOLD=30
        ALT_THRESHOLD=1
        REMOVE_MULTIALLELIC="True"
        PASS="True"
        REMOVE_MITOCHONDRIAL_DNA="True"
        BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
        
        for folder in  ${OUTPUT_VCF_GZ%/*}   ${OUTPUT_FMC_PATH%/*}  ${VEP_FMC_HF_PATH%/*}  ; do
            if [ ! -d $folder ] ; then
                mkdir $folder
            fi
        done
        #01. Mutect2 call
        # qsub -pe smp 6 -e $logPath"/01.MTcall" -o $logPath"/01.MTcall" -N 'MT_01..'${Sample_ID}"_"${TISSUE} -hold_jid "doc_"${Sample_ID}"_"${TISSUE} ${CURRENT_PATH}"/mutect_pipe_01.call.sh" \
        # --Sample_ID ${Sample_ID} --CASE_BAM_PATH ${CASE_BAM_PATH} --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} \
        # --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ}  \
        # --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH}

        
        #02. FMC & HF & RMBLACK
        # qsub -pe smp 5 -e $logPath"/02.FMC_HF_RMBLACK" -o $logPath"/02.FMC_HF_RMBLACK" -N 'MT_02..'${Sample_ID}"_"${TISSUE} -hold_jid  'MT_01..'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/mutect_pipe_02.FMC_HF_RMBLACK.sh" \
        # --Sample_ID ${Sample_ID} --CASE_BAM_PATH ${CASE_BAM_PATH} --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} \
        # --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} --OUTPUT_FMC_PATH ${OUTPUT_FMC_PATH} --OUTPUT_FMC_HF_PATH ${OUTPUT_FMC_HF_PATH}  --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
        # --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH} \
        # --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} \
        # --BLACKLIST ${BLACKLIST}

        # 03. VEP annotation + Nearest gene annotation 
        INPUT_VCF=${OUTPUT_FMC_HF_RMBLACK_PATH}
        VEP_FMC_HF_RMBLACK_PATH=${DATA_PATH%/*}"/04.mutect/03.vep/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf"
        # qsub -pe smp 6 -e $logPath"/03.vep" -o $logPath"/03.vep" -N "MT_03.."${Sample_ID}"_"${TISSUE} -hold_jid 'MT_02..'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/mutect_pipe_20.vep.sh" \
        # --REF ${REF} --INPUT_VCF ${INPUT_VCF} --OUTPUT_VCF ${VEP_FMC_HF_RMBLACK_PATH}      

    done
done




######################################  MULTIPLE CALL ########################################

for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102
    CASE_BAM_PATH1=${DATA_PATH}"/"${hg}"/Tumor/05.Final_bam/"${Sample_ID}"_Tumor.bam"
    CASE_BAM_PATH2=${DATA_PATH}"/"${hg}"/Dura/05.Final_bam/"${Sample_ID}"_Dura.bam"
    CASE_BAM_PATH3=${DATA_PATH}"/"${hg}"/Ventricle/05.Final_bam/"${Sample_ID}"_Ventricle.bam"
    if [ -f ${CASE_BAM_PATH3} ]; then     # File이 있어야만 진행
        CASE_BAM_PATH3="None"
    fi
    CASE_BAM_PATH4=${DATA_PATH}"/"${hg}"/Cortex/05.Final_bam/"${Sample_ID}"_Cortex.bam"
    if [ -f ${CASE_BAM_PATH4} ]; then     # File이 있어야만 진행
        CASE_BAM_PATH4="None"
    fi
    
    CONTROL_BAM_PATH=${DATA_PATH}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
    OUTPUT_VCF_GZ=${DATA_PATH%/*}"/04.mutect/01.raw/"${Sample_ID}"_multiple.vcf.gz"
    OUTPUT_FMC_PATH=${DATA_PATH%/*}"/04.mutect/02.PASS/"${Sample_ID}"_multiple.MT2.FMC.vcf"
    OUTPUT_FMC_HF_PATH=${DATA_PATH%/*}"/04.mutect/02.PASS/"${Sample_ID}"_multiple.MT2.FMC.HF.vcf"
    OUTPUT_FMC_HF_RMBLACK_PATH=${DATA_PATH%/*}"/04.mutect/02.PASS/"${Sample_ID}"_multiple.MT2.FMC.HF.RMBLACK.vcf"
    SAMPLE_THRESHOLD="Dura,Tumor"   # "all"
    DP_THRESHOLD=30
    ALT_THRESHOLD=2
    REMOVE_MULTIALLELIC="True"
    PASS="True"
    REMOVE_MITOCHONDRIAL_DNA="True"
    BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
    
    #11. Multiple sample Mutect2 call (FMC, HF, RMBLACK 포함)
    qsub -pe smp 6 -e $logPath"/11.multipleMT" -o $logPath"/11.multipleMT" -N "MT_11.."${Sample_ID}"_multiple" -hold_jid "doc_"${Sample_ID}"_Dura,doc_"${Sample_ID}"_Tumor" ${CURRENT_PATH}"/mutect_pipe_11.multipleMT.sh" \
    --Sample_ID ${Sample_ID} \
    --CASE_BAM_PATH1 ${CASE_BAM_PATH1} --CASE_BAM_PATH2 ${CASE_BAM_PATH2}  --CASE_BAM_PATH3 ${CASE_BAM_PATH3} --CASE_BAM_PATH4 ${CASE_BAM_PATH4}  \
     --CONTROL_BAM_PATH ${CONTROL_BAM_PATH} \
    --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} --OUTPUT_FMC_PATH ${OUTPUT_FMC_PATH}  --OUTPUT_FMC_HF_PATH ${OUTPUT_FMC_HF_PATH} --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
    --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH} \
    --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} \
    --BLACKLIST ${BLACKLIST}
    
    #12. VEP annotation  +  Nearest gene annotation
    INPUT_VCF=${OUTPUT_FMC_HF_RMBLACK_PATH}
    VEP_FMC_HF_RMBLACK_PATH=${DATA_PATH%/*}"/04.mutect/03.vep/"${Sample_ID}"_multiple.MT2.FMC.HF.RMBLACK.vep.vcf"
    # qsub -pe smp 6 -e $logPath"/12.vep" -o $logPath"/12.vep" -N "MT_12.."${Sample_ID}"_multiple" -hold_jid "MT_11.."${Sample_ID}"_multiple"  ${CURRENT_PATH}"/mutect_pipe_20.vep.sh" \
    # --REF ${REF} --INPUT_VCF ${INPUT_VCF} --OUTPUT_VCF ${VEP_FMC_HF_RMBLACK_PATH}

done




############################################### RESCUE ########################################


# for idx in ${!sample_name_LIST[@]}; do
#     Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102
    
#     for TISSUE in Tumor Dura ; do   #Tumor Dura
#         MULTIPLE_VCF_PATH=${DATA_PATH%/*}"/04.mutect/03.vep/"${Sample_ID}"_multiple.MT2.FMC.HF.RMBLACK.vep.vcf"
#         MULTIPLE_VCF_GZ_PATH=${DATA_PATH%/*}"/04.mutect/03.vep/"${Sample_ID}"_multiple.MT2.FMC.HF.RMBLACK.vep.vcf.gz"
#         INDIVIDUAL_VCF_PATH=${DATA_PATH%/*}"/04.mutect/03.vep/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf"
#         INDIVIDUAL_VCF_GZ_PATH=${DATA_PATH%/*}"/04.mutect/03.vep/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf.gz"
#         INDIVIDUAL_RESCUED_VCF_PATH=${DATA_PATH%/*}"/04.mutect/04.rescue/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.rescue.vcf"
#         INDIVIDUAL_RESCUED_VCF_GZ_PATH=${DATA_PATH%/*}"/04.mutect/04.rescue/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.rescue.vcf.gz"
#         INDIVIDUAL_UNIQUE_VCF_PATH=${DATA_PATH%/*}"/04.mutect/05.unique/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.rescue.unique.vcf"
#         INDIVIDUAL_UNIQUE_VCF_GZ_PATH=${DATA_PATH%/*}"/04.mutect/05.unique/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.rescue.unique.vcf.gz"

#         for dir in ${INDIVIDUAL_RESCUED_VCF_PATH%/*}  ${INDIVIDUAL_UNIQUE_VCF_PATH%/*} ; do
#             if [ ! -d ${dir} ] ; then
#                 mkdir ${dir}
#             fi
#         done

#         # 04.rescue & 05.unique
#         qsub -pe smp 1 -e $logPath"/04.rescue" -o $logPath"/04.rescue" -N "MT_res."${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/mutect_pipe_21.rescue.sh" \
#             --Sample_ID ${Sample_ID} --TISSUE ${TISSUE} \
#             --MULTIPLE_VCF_PATH ${MULTIPLE_VCF_PATH}  --MULTIPLE_VCF_GZ_PATH ${MULTIPLE_VCF_GZ_PATH} \
#             --INDIVIDUAL_VCF_PATH ${INDIVIDUAL_VCF_PATH} --INDIVIDUAL_VCF_GZ_PATH ${INDIVIDUAL_VCF_GZ_PATH} \
#             --INDIVIDUAL_RESCUED_VCF_PATH ${INDIVIDUAL_RESCUED_VCF_PATH} --INDIVIDUAL_RESCUED_VCF_GZ_PATH ${INDIVIDUAL_RESCUED_VCF_GZ_PATH} \
#             --INDIVIDUAL_UNIQUE_VCF_PATH ${INDIVIDUAL_UNIQUE_VCF_PATH} --INDIVIDUAL_UNIQUE_VCF_GZ_PATH ${INDIVIDUAL_UNIQUE_VCF_GZ_PATH}


#         for dir in ${DATA_PATH%/*}"/04.mutect/06.maf/01.shared_yes/" ${DATA_PATH%/*}"/04.mutect/06.maf/02.shared_no/" ; do
#             if [ ! -d ${dir} ] ; then
#                 mkdir -p ${dir}
#             fi
#         done

#         # 06. MAF annotation
#         INPUT_VCF=${INDIVIDUAL_RESCUED_VCF_PATH}
#         OUTPUT_MAF=${DATA_PATH%/*}"/04.mutect/06.maf/01.shared_yes/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.maf"
#         qsub -pe smp 1 -e $logPath"/06.maf" -o $logPath"/06.maf" -N "MT_maf."${Sample_ID}"_"${TISSUE} -hold_jid "MT_res."${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/mutect_pipe_22.manualmaf.sh" \
#         --INPUT_VCF ${INPUT_VCF} --OUTPUT_MAF ${OUTPUT_MAF} --SELECTED_DB "None"

#         INPUT_VCF=${INDIVIDUAL_UNIQUE_VCF_PATH}
#         OUTPUT_MAF=${DATA_PATH%/*}"/04.mutect/06.maf/02.shared_no/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.maf"
#         qsub -pe smp 1 -e $logPath"/06.maf" -o $logPath"/06.maf" -N "MT_maf."${Sample_ID}"_"${TISSUE} -hold_jid "MT_res."${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/mutect_pipe_22.manualmaf.sh" \
#         --INPUT_VCF ${INPUT_VCF} --OUTPUT_MAF ${OUTPUT_MAF} --SELECTED_DB "None"
#     done
# done