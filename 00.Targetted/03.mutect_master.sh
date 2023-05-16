#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/00.Targetted"


PON="/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
#REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
REF="/home/goldpm1/reference/genome.fa"
hg="hg38"
gnomad="/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
TMP_PATH=${DATA_PATH}"/temp"
INTERVAL="/home/goldpm1/resources/TMB359.theragen.hg38.bed"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi

for sublog in 21.MTcall 22.FMC_HF_RMBLACK 23.split 24.MF 25.formatchange 26.vep; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


sample_name_list=$(cat ${CURRENT_PATH}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬


for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102
    
    for TISSUE in Dura ; do   #Tumor Dura
        CASE_BAM_DIR=${DATA_PATH}"/02.Align/"${hg}"/"${TISSUE}"/05.Final_bam"
        CASE_BAM_PATH=${DATA_PATH}"/02.Align/"${hg}"/"${TISSUE}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
        
        OUTPUT_VCF_GZ=${DATA_PATH}"/06.Mutect/01.raw/"${Sample_ID}"_"${TISSUE}".vcf.gz"
        OUTPUT_FMC_PATH=${DATA_PATH}"/06.Mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.vcf"
        OUTPUT_FMC_HF_PATH=${DATA_PATH}"/06.Mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.vcf"
        OUTPUT_FMC_HF_RMBLACK_PATH=${DATA_PATH}"/06.Mutect/02.PASS/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vcf"
        OUTPUT_SNP_VCF=${DATA_PATH}"/06.Mutect/11.split/"${Sample_ID}"_"${TISSUE}".MT2.snp.vcf"
        OUTPUT_INDEL_VCF=${DATA_PATH}"/06.Mutect/11.split/"${Sample_ID}"_"${TISSUE}".MT2.indel.vcf"
        MF_DIR=${DATA_PATH}"/06.Mutect/12.MF"
        VEP_FMC_HF_RMBLACK_PATH=${DATA_PATH}"/06.Mutect/20.vep/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.vcf"

        SAMPLE_THRESHOLD="Dura"   # "all"
        DP_THRESHOLD=30
        ALT_THRESHOLD=1
        REMOVE_MULTIALLELIC="True"
        PASS="True"
        REMOVE_MITOCHONDRIAL_DNA="True"
        BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
        
        for folder in  ${OUTPUT_VCF_GZ%/*}  ${OUTPUT_SNP_VCF%/*}  ${OUTPUT_FMC_PATH%/*}  ${VEP_FMC_HF_RMBLACK_PATH%/*} ${MF_DIR} ; do
            if [ ! -d $folder ] ; then
                mkdir $folder
            fi
        done

        #01. Mutect2 call
        # qsub -pe smp 6 -e $logPath"/21.MTcall" -o $logPath"/21.MTcall" -N 'MT_01.'${Sample_ID}"_"${TISSUE} -hold_jid "doc_"${Sample_ID}"_"${TISSUE} ${CURRENT_PATH}"/03.mutect_pipe_01.call.sh" \
        # --Sample_ID ${Sample_ID} --CASE_BAM_PATH ${CASE_BAM_PATH}  \
        # --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ}  \
        # --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH}

        
        #02. FMC & HF & RMBLACK
        # qsub -pe smp 5 -e $logPath"/22.FMC_HF_RMBLACK" -o $logPath"/22.FMC_HF_RMBLACK" -N 'MT_02.'${Sample_ID}"_"${TISSUE} -hold_jid  'MT_01.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/03.mutect_pipe_02.FMC_HF_RMBLACK.sh" \
        # --Sample_ID ${Sample_ID} --CASE_BAM_PATH ${CASE_BAM_PATH}  \
        # --OUTPUT_VCF_GZ ${OUTPUT_VCF_GZ} --OUTPUT_FMC_PATH ${OUTPUT_FMC_PATH} --OUTPUT_FMC_HF_PATH ${OUTPUT_FMC_HF_PATH}  --OUTPUT_FMC_HF_RMBLACK_PATH ${OUTPUT_FMC_HF_RMBLACK_PATH} \
        # --PON ${PON} --REF ${REF} --gnomad ${gnomad} --INTERVAL ${INTERVAL} --TMP_PATH ${TMP_PATH} \
        # --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD} --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA} \
        # --BLACKLIST ${BLACKLIST}

        #11. snp, indel split (for MF)
        INPUT_VCF_GZ=${OUTPUT_VCF_GZ}
        OUTPUT_SNP_VCF=${DATA_PATH}"/06.Mutect/11.split/"${Sample_ID}"_"${TISSUE}".MT2.snp.vcf"
        OUTPUT_INDEL_VCF=${DATA_PATH}"/06.Mutect/11.split/"${Sample_ID}"_"${TISSUE}".MT2.indel.vcf"
        # qsub -pe smp 2 -e $logPath"/23.split" -o $logPath"/23.split" -N 'MT_11.'${Sample_ID}"_"${TISSUE} -hold_jid "MT_01."${Sample_ID}"_"${TISSUE} ${CURRENT_PATH}"/03.mutect_pipe_11.snpindelsplit.sh" \
        #     --INPUT_VCF_GZ ${INPUT_VCF_GZ} \
        #     --OUTPUT_SNP_VCF ${OUTPUT_SNP_VCF} \
        #     --OUTPUT_INDEL_VCF ${OUTPUT_INDEL_VCF}
        
        #12. MF 
        MT_SNP_VCF=${OUTPUT_SNP_VCF}
        MT_INDEL_VCF=${OUTPUT_INDEL_VCF}
        qsub -pe smp 2 -e $logPath"/24.MF" -o $logPath"/24.MF" -N 'MT_12.'${Sample_ID}"_"${TISSUE} -hold_jid "MT_11."${Sample_ID}"_"${TISSUE} ${CURRENT_PATH}"/03.mutect_pipe_12.MF.sh" \
            --REF ${REF} \
            --BAM_DIR ${CASE_BAM_DIR} \
            --CASE_BAM_PATH ${CASE_BAM_PATH} \
            --ID ${Sample_ID}"_"${TISSUE}  \
            --MT_SNP_VCF ${MT_SNP_VCF} \
            --MT_INDEL_VCF ${MT_INDEL_VCF} \
            --MF_DIR ${MF_DIR} \
            --F_bed "/data/project/MRS/Resource/MF/fout_snv.bed" \
            --F_bed_ind "/data/project/MRS/Resource/MF/fout_ind.bed" \
            --MF_TOOL_DIR "/home/goldpm1/tools/MosaicForecast2"



        #14. Bed 파일로 만들고 Mutect2 vcf와 bedtools intersect하기 (원래 폴더에 final.vcf.gz 생성)
        # AnalysisPath="/data/project/craniosynostosis/9.mosaic/mf"
        # MT_Path="/data/project/craniosynostosis/9.mosaic/mutect2"
        # MFSNP=$AnalysisPath"/"$ID"/"$ID".SNP.Predictions.mosaic"
        # MFIND=$AnalysisPath"/ind."$ID"/"$ID".IND.Predictions.mosaic"

        # 20. VEP annotation + Nearest gene annotation   (이거는 base로 하는 게 좋다. conda activate cnvpytor 하면 가끔 에러난다)
        # INPUT_VCF=${OUTPUT_FMC_HF_RMBLACK_PATH}
        # qsub -pe smp 6 -e $logPath"/26.vep" -o $logPath"/26.vep" -N "MT_03."${Sample_ID}"_"${TISSUE} -hold_jid 'MT_02.'${Sample_ID}"_"${TISSUE}  ${CURRENT_PATH}"/mutect_pipe_20.vep.sh" \
        # --REF ${REF} --INPUT_VCF ${INPUT_VCF} --OUTPUT_VCF ${VEP_FMC_HF_RMBLACK_PATH}      

    done
done
    
