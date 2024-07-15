#!/bin/bash
#$ -cwd
#$ -S /bin/bash

CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/home/goldpm1/Meningioma/02.Align"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "mutTR_01.bySequenza" "mutTR_02.byFacetCNV"; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done

REF="/home/goldpm1/reference/genome.fa"
#REF="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
hg="hg38"
INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz"
    
MUTATIONTIMER_DIR="/data/project/Meningioma/31.Clonality/03.mutationtimeR"
MUTECT_VCF_DIR="/data/project/Meningioma/04.mutect/04.rescue/"
SEQUENZA_DIR="/data/project/Meningioma/11.cnv/2.sequenza"
FACET_CNV_DIR="/data/project/Meningioma/31.Clonality/01.make_matrix"

MUTATIONTIMER_INPUT_VCF_DIR=${MUTATIONTIMER_DIR}"/01.vcf"

for dir in $MUTATIONTIMER_INPUT_VCF_DIR $MUTATIONTIMER_RESULT_DIR; do
    if [ ! -d ${dir} ] ; then
        mkdir -p ${dir}
    fi
done


sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
sample_name_LIST=(${sample_name_list// / })     # array로 만듬



for idx in ${!sample_name_LIST[@]}; do
    Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102
    for TISSUE in Dura Tumor; do
        echo -e ${Sample_ID}"_"${TISSUE}

        #01.  Mutect call -> Blood를 samplename에서 빼기  +   chr을 빼서  한 곳으로 옮기기
        MUTECT_VCF_PATH=${MUTECT_VCF_DIR}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.rescue.vcf"
        MUTATIONTIMER_INPUT_VCF_PATH=${MUTATIONTIMER_INPUT_VCF_DIR}"/"${Sample_ID}"_"${TISSUE}".MT2.FMC.HF.RMBLACK.vep.rescue.removechr.vcf"

        # bcftools view -s ${Sample_ID}"_"${TISSUE} ${MUTECT_VCF_PATH} > ${MUTATIONTIMER_INPUT_VCF_PATH}".temp"  # Blood 를 vcf에서 뺴기
        # grep -v 'chr.*random' ${MUTATIONTIMER_INPUT_VCF_PATH}".temp" | grep -v 'chr.*_alt' |  grep -v "chrM" | grep -v "chrEBV"| grep -v "chrUn" > ${MUTATIONTIMER_INPUT_VCF_PATH}".temp2"
        # awk '{gsub(/^chr/,""); print}' ${MUTATIONTIMER_INPUT_VCF_PATH}".temp2" > ${MUTATIONTIMER_INPUT_VCF_PATH}".temp3"  # variant에서 chr1 -> 1 로 변형해주기
        # sed -i 's/ID=chr/ID=/g' ${MUTATIONTIMER_INPUT_VCF_PATH}".temp3"  # contig에서 ID=chr1 -> ID= 로 변형해주기
        # mv ${MUTATIONTIMER_INPUT_VCF_PATH}".temp3" ${MUTATIONTIMER_INPUT_VCF_PATH}
        # rm -rf ${MUTATIONTIMER_INPUT_VCF_PATH}".temp" ${MUTATIONTIMER_INPUT_VCF_PATH}".temp2" ${MUTATIONTIMER_INPUT_VCF_PATH}".temp3" ${MUTATIONTIMER_INPUT_VCF_PATH}".temp4"
    
        
        #02. MUTATIONTIMER 실행하기 (Sequenza)
        SEQUENZA_SEGMENT_PATH=${SEQUENZA_DIR}"/hg19to38/"${Sample_ID}"_"${TISSUE}"_segments.txt"
        SEQUENZA_PLOIDY_PATH=${SEQUENZA_DIR}"/hg19/"${Sample_ID}"_"${TISSUE}"_confints_CP.txt"             # 이건 hg19에만 있다
        MUTATIONTIMER_RESULT_DIR=${MUTATIONTIMER_DIR}"/02.result_sequenza"
        if [ ! -d ${MUTATIONTIMER_RESULT_DIR} ] ; then
            mkdir -p ${MUTATIONTIMER_RESULT_DIR}
        fi
        qsub -pe smp 2 -e $logPath"/mutTR_01.bySequenza" -o $logPath"/mutTR_01.bySequenza" -N "mutTR_"${Sample_ID}"_"${TISSUE}"_bySequenza"  ${CURRENT_PATH}"/MutationTimeR_pipe_01.bySequenza.sh" \
            --Sample_ID ${Sample_ID} \
            --TISSUE ${TISSUE} \
            --SEQUENZA_SEGMENT_PATH ${SEQUENZA_SEGMENT_PATH} \
            --SEQUENZA_PLOIDY_PATH ${SEQUENZA_PLOIDY_PATH} \
            --MUTATIONTIMER_INPUT_VCF_PATH ${MUTATIONTIMER_INPUT_VCF_PATH} \
            --MUTATIONTIMER_RESULT_DIR ${MUTATIONTIMER_RESULT_DIR}

        #03. MUTATIONTIMER 실행하기 (FacetCNV)
        FACET_CNV_MATRIX_PATH="/data/project/Meningioma/31.Clonality/01.make_matrix/"${Sample_ID}"/"${Sample_ID}"_"${TISSUE}".facetcnv_to_bed_df.tsv"
        MUTATIONTIMER_RESULT_DIR=${MUTATIONTIMER_DIR}"/02.result_facetcnv"
        if [ ! -d ${MUTATIONTIMER_RESULT_DIR} ] ; then
            mkdir -p ${MUTATIONTIMER_RESULT_DIR}
        fi
        qsub -pe smp 2 -e $logPath"/mutTR_02.byFacetCNV" -o $logPath"/mutTR_02.byFacetCNV" -N "mutTR_"${Sample_ID}"_"${TISSUE}"_byFacetCNV"  ${CURRENT_PATH}"/MutationTimeR_pipe_02.byFacetCNV.sh" \
            --Sample_ID ${Sample_ID} \
            --TISSUE ${TISSUE} \
            --FACET_CNV_MATRIX_PATH ${FACET_CNV_MATRIX_PATH} \
            --MUTATIONTIMER_INPUT_VCF_PATH ${MUTATIONTIMER_INPUT_VCF_PATH} \
            --MUTATIONTIMER_RESULT_DIR ${MUTATIONTIMER_RESULT_DIR}
    done
done
