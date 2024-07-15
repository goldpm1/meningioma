#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#PROJECT_DIR="/data/project/Meningioma/61.Lowinput/01.XT_HS"
PROJECT_DIR="/data/project/Meningioma/61.Lowinput/02.PTA"

WES_BAM_DIR="/data/project/Meningioma/02.Align/"
BAM_DIR=${PROJECT_DIR}"/02.Align"
MUTECT_DIR="/data/project/Meningioma/04.mutect"
GVCF_DIR=${PROJECT_DIR}"/05.gvcf"
HC_DIR=${PROJECT_DIR}"/06.HC"
SCRIPT_DIR="/data/project/Meningioma/script/61.Lowinput"
logPath=${SCRIPT_DIR}"/log"
PON="/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
#REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
REF="/home/goldpm1/reference/genome.fa"
hg="hg38"
gnomad="/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz" 
BLACKLIST="/home/goldpm1/resources/RM+SegDup.bed"
TMP_PATH=${BAM_DIR}"/temp"
INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
INTERVAL="/home/goldpm1/resources/whole.exome.bed"




for sublog in "06-1.call" "06-2.VariantFiltration" "06-3.combine_genotype_gvcf" "06-4.VQSR" "06-5.pick_somatic"; do
    if [ -d ${logPath}"/"${sublog} ]; then
        rm -rf ${logPath}"/"${sublog}
    fi
    if [ ! -d ${logPath}"/"${sublog} ]; then
        mkdir -p ${logPath}"/"${sublog}
    fi
done
for subdir in "01.call" "02.VariantFiltration" "03.combine_genotype" "04.VQSR" "05.pick_somatic"  ; do
    if [ ! -d ${HC_DIR}"/"${subdir} ]; then
        mkdir -p ${HC_DIR}"/"${subdir}
    fi
done



# 01. call (single)
#for Sample_ID in 190426; do
for Sample_ID in 230405; do
    HC_VCF_LIST=""
    BAM_DIR_LIST=""
    TITLE_LIST=""
    hold_j=""

    for Clone_No in 1 2 3 4 5 6 7 8 ; do
        CASE_BAM_PATH=${BAM_DIR}"/"${hg}"/Tumor/06.Final_bam/"${Sample_ID}"_"${Clone_No}".bam"
        BAM_DIR_LIST=${BAM_DIR_LIST}","${CASE_BAM_PATH}
        TITLE_LIST=${TITLE_LIST}","${Sample_ID}"_"${Clone_No}
        echo -e ${Sample_ID}"_"${Clone_No}
            
        # # [HC call]
        # qsub -pe smp 5 -e $logPath"/06-1.call" -o $logPath"/06-1.call"  -hold_jid "postbwa_"${Sample_ID}"_"${Clone_No} -N 'HC1_'${Sample_ID}"_"${Clone_No} ${SCRIPT_DIR}"/06.HC_pipe_01.call.sh" \
        #         --BAM_PATH ${CASE_BAM_PATH} \
        #         --INTERVAL ${INTERVAL} \
        #         --REF ${REF} \
        #         --dbSNP ${dbSNP} \
        #         --OUTPUT_HC_GVCF ${HC_DIR}"/01.call/"${Sample_ID}"_"${Clone_No}".hc.gvcf" \
        #         --OUTPUT_HC_VCF ${HC_DIR}"/01.call/"${Sample_ID}"_"${Clone_No}".hc.vcf"
                
        # #[VariantFiltration by INFO filed]  
        # qsub -pe smp 5 -e $logPath"/06-2.VariantFiltration" -o $logPath"/06-2.VariantFiltration" -hold_jid "HC1_"${Sample_ID}"_"${Clone_No} -N 'HC2_gvcf_'${Sample_ID}"_"${Clone_No} ${SCRIPT_DIR}"/06.HC_pipe_02.VariantFiltration.sh" \
        #     --REF ${REF} \
        #     --BLACKLIST ${BLACKLIST} \
        #     --dbSNP "False" \
        #     --OUTPUT_HC ${HC_DIR}"/01.call/"${Sample_ID}"_"${Clone_No}".hc.gvcf" \
        #     --OUTPUT_HC_VF ${HC_DIR}"/02.VariantFiltration/"${Sample_ID}"_"${Clone_No}".hc.VF.gvcf" 

        # qsub -pe smp 5 -e $logPath"/06-2.VariantFiltration" -o $logPath"/06-2.VariantFiltration" -hold_jid "HC1_"${Sample_ID}"_"${Clone_No} -N 'HC2_vcf_'${Sample_ID}"_"${Clone_No} ${SCRIPT_DIR}"/06.HC_pipe_02.VariantFiltration.sh" \
        #     --REF ${REF} \
        #     --BLACKLIST ${BLACKLIST} \
        #     --dbSNP ${dbSNP} \
        #     --OUTPUT_HC ${HC_DIR}"/01.call/"${Sample_ID}"_"${Clone_No}".hc.vcf" \
        #     --OUTPUT_HC_VF ${HC_DIR}"/02.VariantFiltration/"${Sample_ID}"_"${Clone_No}".hc.VF.vcf"


        hold_j=${hold_j}",HC2_gvcf_"${Sample_ID}"_"${Clone_No}
        HC_VCF_LIST=${HC_VCF_LIST}","${HC_DIR}"/02.VariantFiltration/"${Sample_ID}"_"${Clone_No}".hc.VF.gvcf.gz"
    done
    
    
    ############################ Blood도 GVCF 해주자 ############################
    WES_TUMOR_BAM_PATH=${WES_BAM_DIR}"/"${hg}"/Tumor/05.Final_bam/"${Sample_ID}"_Tumor.bam"
    if [ ${Sample_ID} == "190426" ]; then
        WES_TUMOR_BAM_PATH=${WES_BAM_DIR}"/"${hg}"/Tumor/05.Final_bam/"${Sample_ID}"_Tumor_PT.bam"
    fi
    WES_TUMOR_BED="/data/project/Meningioma/04.mutect/03.Tumor_interval/"${Sample_ID}"_Tumor.MT2.FMC.HF.RMBLACK.bed"
    WES_DURA_BAM_PATH=${WES_BAM_DIR}"/"${hg}"/Dura/05.Final_bam/"${Sample_ID}"_Dura.bam"
    CONTROL_BAM_PATH=${WES_BAM_DIR}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
    BAM_DIR_LIST=${BAM_DIR_LIST}","${WES_TUMOR_BAM_PATH}","${WES_DURA_BAM_PATH}","${CONTROL_BAM_PATH}
    TITLE_LIST=${TITLE_LIST}","${Sample_ID}"_Tumor,"${Sample_ID}"_Dura,"${Sample_ID}"_Blood"
    echo -e ${Sample_ID}"_Blood"
    # # [HC call]
    # qsub -pe smp 5 -e $logPath"/06-1.call" -o $logPath"/06-1.call"  -N 'HC1_'${Sample_ID}"_Blood" ${SCRIPT_DIR}"/06.HC_pipe_01.call.sh" \
    #         --BAM_PATH ${CONTROL_BAM_PATH} \
    #         --INTERVAL ${INTERVAL} \
    #         --REF ${REF} \
    #         --dbSNP ${dbSNP} \
    #         --OUTPUT_HC_GVCF ${HC_DIR}"/01.call/"${Sample_ID}"_Blood.hc.gvcf" \
    #         --OUTPUT_HC_VCF ${HC_DIR}"/01.call/"${Sample_ID}"_Blood.hc.vcf"
            
    # #[VariantFiltration by INFO filed]  
    # qsub -pe smp 5 -e $logPath"/06-2.VariantFiltration" -o $logPath"/06-2.VariantFiltration" -hold_jid "HC1_"${Sample_ID}"_Blood" -N 'HC2_gvcf_'${Sample_ID}"_Blood" ${SCRIPT_DIR}"/06.HC_pipe_02.VariantFiltration.sh" \
    #     --REF ${REF} \
    #     --BLACKLIST ${BLACKLIST} \
    #     --dbSNP "False" \
    #     --OUTPUT_HC ${HC_DIR}"/01.call/"${Sample_ID}"_Blood.hc.gvcf" \
    #     --OUTPUT_HC_VF ${HC_DIR}"/02.VariantFiltration/"${Sample_ID}"_Blood.hc.VF.gvcf"

    # qsub -pe smp 5 -e $logPath"/06-2.VariantFiltration" -o $logPath"/06-2.VariantFiltration" -hold_jid "HC1_"${Sample_ID}"_Blood" -N 'HC2_vcf_'${Sample_ID}"_Blood" ${SCRIPT_DIR}"/06.HC_pipe_02.VariantFiltration.sh" \
    #     --REF ${REF} \
    #     --BLACKLIST ${BLACKLIST} \
    #     --dbSNP ${dbSNP} \
    #     --OUTPUT_HC ${HC_DIR}"/01.call/"${Sample_ID}"_Blood.hc.vcf" \
    #     --OUTPUT_HC_VF ${HC_DIR}"/02.VariantFiltration/"${Sample_ID}"_Blood.hc.VF.vcf" 


    hold_j=${hold_j}",HC2_"${Sample_ID}"_Blood"
    HC_VCF_LIST=${HC_VCF_LIST}","${HC_DIR}"/02.VariantFiltration/"${Sample_ID}"_Blood.hc.VF.gvcf.gz"


    #03. Combine/GenomicsDBImport & Genotype GVCF (위에꺼 다 끝나고 진행; conda deactivate)
    COMBINED_GVCF=${HC_DIR}"/03.combine_genotype/"${Sample_ID}".combined.vcf.gz"
    GENOMICSDB_WORKSPACE=${HC_DIR}"/03.combine_genotype/"${Sample_ID}".genomicsDB"
    GENOTYPE_GVCF=${HC_DIR}"/03.combine_genotype/"${Sample_ID}".genotype.vcf.gz"
    INTERVAL="/home/goldpm1/resources/genomicDBImport/S04380110_Covered.100pad.merged.bed"
    #INTERVAL=""/home/goldpm1/resources/whole.chromosome.GRCh38.bed""
    
    # python3 ${SCRIPT_DIR}"/06.HC_pipe_03.combine_genotype_gvcf.py" \
    #     --REF ${REF} \
    #     --dbSNP ${dbSNP} \
    #     --HC_VCF_LIST "${HC_VCF_LIST#?}" \
    #     --SAMPLE_ID ${Sample_ID} \
    #     --SCRIPT_DIR ${SCRIPT_DIR} \
    #     --COMBINED_GVCF ${COMBINED_GVCF} \
    #     --GENOMICSDB_WORKSPACE ${GENOMICSDB_WORKSPACE} \
    #     --GENOTYPE_GVCF ${GENOTYPE_GVCF} \
    #     --LOGPATH ${logPath}"/06-3.combine_genotype_gvcf" \
    #     --TMP_DIR ${TMP_PATH} \
    #     --INTERVAL ${INTERVAL}


    #04. VQSR
    qsub -pe smp 5 -e $logPath"/06-4.VQSR" -o $logPath"/06-4.VQSR" -hold_jid HC3_${Sample_ID} -N 'HC4_'${Sample_ID} ${SCRIPT_DIR}"/06.HC_pipe_04.VQSR.sh" \
        --REF ${REF} \
        --INPUT_VCF ${GENOTYPE_GVCF} \
        --INPUT_SNP_VCF ${GENOTYPE_GVCF%vcf.gz}"snp.vcf.gz" \
        --INPUT_INDEL_VCF ${GENOTYPE_GVCF%vcf.gz}"indel.vcf.gz" \
        --OUTPUT_SNP_VCF ${HC_DIR}"/04.VQSR/"${Sample_ID}".VQSR.snp.vcf" \
        --OUTPUT_INDEL_VCF ${HC_DIR}"/04.VQSR/"${Sample_ID}".VQSR.indel.vcf"  \
        --OUTPUT_VCF ${HC_DIR}"/04.VQSR/"${Sample_ID}".VQSR.vcf"  \
        --TRANCHES_SNP ${HC_DIR}"/04.VQSR/"${Sample_ID}".snp.tranches"  \
        --RECAL_SNP ${HC_DIR}"/04.VQSR/"${Sample_ID}".snp.recal"  \
        --TRANCHES_INDEL ${HC_DIR}"/04.VQSR/"${Sample_ID}".indel.tranches"  \
        --RECAL_INDEL ${HC_DIR}"/04.VQSR/"${Sample_ID}".indel.tranches"

    VQSR_VCF_PATH=${HC_DIR}"/04.VQSR/"${Sample_ID}".VQSR.vcf"

    #05. Pick shared germline variants
    OUTPUT_VCF=${HC_DIR}"/05.pick_somatic/"${Sample_ID}".pick_somatic.vcf"

    BCFTOOLS_MERGE_TXT=${PROJECT_DIR}"/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/"${Sample_ID}".txt"
    BCFTOOLS_MERGE_OUTPUT_VCF_GZ=${PROJECT_DIR}"/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/"${Sample_ID}".merge.vcf.gz"
    if [ ! -d ${BCFTOOLS_MERGE_OUTPUT_VCF_GZ%/*} ]; then
        mkdir -p ${BCFTOOLS_MERGE_OUTPUT_VCF_GZ%/*}
    fi
    WES_DURA_VCF_PATH=${MUTECT_DIR}"/04.Other_rescue/"${Sample_ID}"_Dura.MT2.FMC.HF.RMBLACK.rescue.vcf.gz"
    WES_DURA_VCF_GZ_PATH=${MUTECT_DIR}"/04.Other_rescue/"${Sample_ID}"_Dura.MT2.FMC.HF.RMBLACK.rescue.vcf.gz"
    THRESHOLD_DEPTH=40
    if [ ${Sample_ID} == "230405" ]; then
        WES_DURA_VCF_PATH=${MUTECT_DIR}"/04.Other_rescue/"${Sample_ID}"_2_Dura.MT2.FMC.HF.RMBLACK.rescue.vcf.gz"
        WES_DURA_VCF_GZ_PATH=${MUTECT_DIR}"/04.Other_rescue/"${Sample_ID}"_2_Dura.MT2.FMC.HF.RMBLACK.rescue.vcf.gz"
        THRESHOLD_DEPTH=15
    fi
    echo "$WES_DURA_VCF_GZ_PATH" > ${BCFTOOLS_MERGE_TXT}
    echo ${OUTPUT_VCF}".gz" >> ${BCFTOOLS_MERGE_TXT}

    OUTPUT_HEATMAP_PATH=${PROJECT_DIR}"/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/"${Sample_ID}".heatmap.pdf"

    qsub -pe smp 1 -e $logPath"/06-5.pick_somatic" -o $logPath"/06-5.pick_somatic" -hold_jid "HC4_"${Sample_ID} -N 'HC5_'${Sample_ID} ${SCRIPT_DIR}"/06.HC_pipe_05.pick_somatic.sh" \
        --INPUT_VCF ${VQSR_VCF_PATH} \
        --OUTPUT_VCF ${OUTPUT_VCF} \
        --OUTPUT_BAMSNAP_DIR ${HC_DIR}"/05.pick_somatic"  \
        --SAMPLE_ID ${Sample_ID} \
        --BAM_DIR_LIST "${BAM_DIR_LIST#?}" \
        --TITLE_LIST "${TITLE_LIST#?}" \
        --WES_TUMOR_BED  ${WES_TUMOR_BED} \
        --WES_DURA_VCF ${WES_DURA_VCF_PATH} \
        --BCFTOOLS_MERGE_TXT ${BCFTOOLS_MERGE_TXT} \
        --BCFTOOLS_MERGE_OUTPUT_VCF_GZ ${BCFTOOLS_MERGE_OUTPUT_VCF_GZ} \
        --OUTPUT_HEATMAP_PATH ${OUTPUT_HEATMAP_PATH} \
        --THRESHOLD_DEPTH ${THRESHOLD_DEPTH}
done



    #${GENOTYPE_GVCF%.gz}
