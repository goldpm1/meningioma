#!/bin/bash
#$ -S /bin/bash
#$ -cwd

REF_hg19="/home/goldpm1/reference/Broadhg19/Homo_sapiens_assembly19.fasta"             # 왜그런지는 모르겠지만 이 reference를 써야 돌아간다
REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
REF_hg38="/home/goldpm1/reference/genome.fa"
#REF_hg38="/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"
REF=${REF_hg19}
hg="hg19"
WIGGZ="/home/goldpm1/resources/sequenza/reference.hg19.genome.fa.gc50Base.wig.gz"


CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/data/project/Meningioma/02.Align"
SEQUENZApath="/data/project/Meningioma/11.cnv"

if [ ! -d $logPath ] ; then
    mkdir $logPath
fi
for sublog in "seq_01.preprocess" "seq_02.runningR" "seq_03.liftover"; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


# sample_name_list=$(cat ${CURRENT_PATH%/*}"/sample_name.txt")
# sample_name_LIST=(${sample_name_list// / })     # array로 만듬

# for idx in ${!sample_name_LIST[@]}; do
#     Sample_ID=${sample_name_LIST[idx]}        #220930, 221026, 221102



#for Sample_ID in 220930 221026 221102 221202 230127 230323_2 230405_2 230419 230526 230822 230920; do
for Sample_ID in 190426; do
    for TISSUE in Tumor Tumor_FFT Tumor_PCT Tumor_PFT Tumor_PP Tumor_PT Dura Ventricle Cortex; do   #Tumor Dura
    #for TISSUE in Dura; do
        CASE_BAM_PATH=${DATA_PATH}"/"${hg}"/"${TISSUE%_*}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"

        if [ -f ${CASE_BAM_PATH} ]; then     # File이 있어야만 진행
            echo $Sample_ID"_"$TISSUE
        
            CONTROL_BAM_PATH=${DATA_PATH}"/"${hg}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"

            SEQUENZA_SEQZ=${SEQUENZApath}"/1.seqz/"${hg}"/1.original/"${Sample_ID}"_"${TISSUE}".seqz.gz"
            SEQUENZA_SMALL_SEQZ=${SEQUENZApath}"/1.seqz/"${hg}"/2.small/"${Sample_ID}"_"${TISSUE}".small.seqz.gz"
                
            for subpath in ${SEQUENZA_SEQZ} ${SEQUENZA_SMALL_SEQZ}; do
                if [ ! -d ${subpath%/*} ]; then
                    mkdir -p ${subpath%/*}
                fi
            done
        
            #"/home/goldpm1/resources/sequenza/reference.hg19.genome.fa.gc50Base.wig.gz"

            #1. preprocess
            # qsub -pe smp 4 -o $logPath"/seq_01.preprocess" -e $logPath"/seq_01.preprocess" -N "sq1_"$hg"_"$Sample_ID"_"$TISSUE -hold_jid "doc_"${Sample_ID}"_"${TISSUE} "sequenza_pipe_01.preprocess.sh"  \
            #     --CASE_BAM ${CASE_BAM_PATH} \
            #     --CONTROL_BAM ${CONTROL_BAM_PATH} \
            #     --REF ${REF}  \
            #     --WIGGZ ${WIGGZ} \
            #     --SEQUENZA_SEQZ ${SEQUENZA_SEQZ}  \
            #     --SEQUENZA_SMALL_SEQZ ${SEQUENZA_SMALL_SEQZ} 
            

            # #2. runningR 
            SEQUENZA_OUTPUT_DIR=${SEQUENZApath}"/2.sequenza/"${hg}
            for subpath in ${SEQUENZA_OUTPUT_DIR}; do
                if [ ! -d ${subpath} ]; then
                    mkdir -p ${subpath}
                fi
            done
            # qsub -pe smp 5 -o $logPath"/seq_02.runningR" -e $logPath"/seq_02.runningR" -N "sq2_"$hg"_"$Sample_ID"_"$TISSUE -hold_jid "sq1_"$hg"_"$Sample_ID"_"$TISSUE "sequenza_pipe_02.runningR.sh" \
            #     --ID $Sample_ID"_"$TISSUE \
            #     --hg $hg \
            #     --SEQUENZA_SMALL_SEQZ ${SEQUENZA_SMALL_SEQZ}".removechrM.gz"  \
            #     --SEQUENZA_OUTPUT_DIR ${SEQUENZA_OUTPUT_DIR}


            #3. hg19 to hg38 liftover
            SEQUENZA_LIFTOVER_OUTPUT_DIR=${SEQUENZApath}"/2.sequenza/hg19to38"
            for subpath in ${SEQUENZA_LIFTOVER_OUTPUT_DIR}; do
                if [ ! -d ${subpath} ]; then
                    mkdir -p ${subpath}
                fi
            done
            qsub -pe smp 1 -o $logPath"/seq_03.liftover" -e $logPath"/seq_03.liftover" -N "sq3_"$hg"_"$Sample_ID"_"$TISSUE -hold_jid "sq2_"$hg"_"$Sample_ID"_"$TISSUE "sequenza_pipe_03.liftover.sh" \
                --ID $Sample_ID"_"$TISSUE \
                --SEQUENZA_OUTPUT_DIR ${SEQUENZA_OUTPUT_DIR} \
                --SEQUENZA_LIFTOVER_OUTPUT_DIR ${SEQUENZA_LIFTOVER_OUTPUT_DIR}
       
        fi
    done
done