#!/bin/bash
#$ -cwd
#$ -S /bin/bash


REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"
REF_hg19="/home/goldpm1/reference/Broadhg19/Homo_sapiens_assembly19.fasta"             # 왜그런지는 모르겠지만 이 reference를 써야 돌아간다
REF_hg38="/home/goldpm1/reference/genome.fa"
hg="hg38"

INTERVAL="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"
dbSNP="/data/public/dbSNP/b154/GRCh38/GCF_000001405.38.re.common.vcf.gz"


CURRENT_PATH=`pwd -P`
logPath=$CURRENT_PATH"/log"

DATA_PATH="/home/goldpm1/Meningioma/02.Align"
HATCHET_PATH="/home/goldpm1/Meningioma/11.cnv/3.hatchet"

###### conda activate cnvpytor 필요


for sublog in "hatchet_01.genotype-snps" "hatchet_02.count-alleles" "hatchet_03.count-reads" "hatchet_04.phase-snps" "hatchet_05.combine-counts" "hatchet_06.cluster-bins" "hatchet_07.plot-bins" "hatchet_08.compute-cn" "hatchet_09.plot-cn"; do
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

    CONTROL_BAM_PATH=${DATA_PATH}"/"${hg38}"/Blood/05.Final_bam/"${Sample_ID}"_Blood.bam"
    CASE_BAM_PATH=""
    SAMPLES=""
    ALLNAMES=${Sample_ID}"_Blood"

    for TISSUE in Tumor Dura ; do   #Tumor Dura
        SAMPLES=${SAMPLES}","${Sample_ID}"_"${TISSUE}
        ALLNAMES=${ALLNAMES}","${Sample_ID}"_"${TISSUE}
        CASE_BAM_PATH=${CASE_BAM_PATH}","${DATA_PATH}"/"${hg38}"/"${TISSUE}"/05.Final_bam/"${Sample_ID}"_"${TISSUE}".bam"
    done

    # make subpath
    HATCHET_PATH_01=${HATCHET_PATH}"/01.genotype-snps"
    HATCHET_PATH_02=${HATCHET_PATH}"/02.count-alleles"
    HATCHET_PATH_03=${HATCHET_PATH}"/03.count-reads"
    HATCHET_PATH_04=${HATCHET_PATH}"/04.phase-snps"
    HATCHET_PATH_05=${HATCHET_PATH}"/05.combine-counts"
    HATCHET_PATH_06=${HATCHET_PATH}"/06.cluster-bins"
    HATCHET_PATH_07=${HATCHET_PATH}"/07.plot-bins"
    HATCHET_PATH_08=${HATCHET_PATH}"/08.compute-cn"
    HATCHET_PATH_09=${HATCHET_PATH}"/09.plot-cn"
    
    for subpath in ${HATCHET_PATH_01} ${HATCHET_PATH_02} ${HATCHET_PATH_03} ${HATCHET_PATH_04} ${HATCHET_PATH_05} ${HATCHET_PATH_06} ${HATCHET_PATH_07} ${HATCHET_PATH_08} ${HATCHET_PATH_09}; do
        if [ ! -d ${subpath} ]; then
            mkdir -p ${subpath}
        fi
    done



    SAMTOOLS="/opt/Yonsei/samtools/1.7/"
    BCFTOOLS="/opt/Yonsei/bcftools/1.7/"
    BGZIP="~/miniconda3/envs/cnvpytor/bin/"
    SHAPEIT="~/miniconda3/envs/cnvpytor/bin/"
    PICARD="/opt/Yonsei/Picard/2.26.4/"
    MINCOV="30"
    MAXCOV="1000"
    PROCESSES="3"
    READQUALITY="20"
    BASEQUALITY="20"


    # 01. genotype-snps 
    OUTPUTSNPS=${HATCHET_PATH_01}"/"${Sample_ID}
    if [ ! -d ${OUTPUTSNPS} ]; then
        mkdir -p ${OUTPUTSNPS}
    fi
    for sublog in "hatchet_01.genotype-snps"; do
        if [ $logPath"/"$sublog ] ; then
            rm -rf $logPath"/"$sublog
        fi
        if [ ! -d $logPath"/"$sublog ] ; then
            mkdir -p $logPath"/"$sublog
        fi
    done
    # qsub -pe smp 3 -o $logPath"/hatchet_01.genotype-snps" -e $logPath"/hatchet_01.genotype-snps" -N "hat01_"$Sample_ID  "hatchet_pipe_01.genotype-snps.sh" \
    # --NORMAL ${CONTROL_BAM_PATH} --REF ${REF_hg38} --OUTPUTSNPS ${OUTPUTSNPS} \
    # --SAMTOOLS ${SAMTOOLS} --BCFTOOLS ${BCFTOOLS} --MINCOV ${MINCOV} --MAXCOV ${MAXCOV} --PROCESSES ${PROCESSES} --READQUALITY ${READQUALITY} --BASEQUALITY ${BASEQUALITY}

    #02. count-alleles (OUTPUTSNPS에 아무것도 없다. 이상하다. WES에서 약 1시간 30분 걸림)
    # SNPS=${OUTPUTSNPS}
    # OUTPUTNORMAL=${HATCHET_PATH_02}"/"${Sample_ID}"/"${Sample_ID}".outputnormal.tsv"
    # OUTPUTTUMORS=${HATCHET_PATH_02}"/"${Sample_ID}"/"${Sample_ID}".outputtumors.tsv"
    # OUTPUTSNPS=${HATCHET_PATH_02}"/"${Sample_ID}"/"      # 이건 있으나 마나한데 왜 쓰라고 하는지...
    # if [ ! -d ${OUTPUTSNPS%"/*"} ]; then
    #     mkdir -p ${OUTPUTSNPS%"/*"}
    # fi
    # for sublog in "hatchet_02.count-alleles"; do
    #     if [ $logPath"/"$sublog ] ; then
    #         rm -rf $logPath"/"$sublog
    #     fi
    #     if [ ! -d $logPath"/"$sublog ] ; then
    #         mkdir -p $logPath"/"$sublog
    #     fi
    # done
    # qsub -pe smp 3 -o $logPath"/hatchet_02.count-alleles" -e $logPath"/hatchet_02.count-alleles" -N "hat02_"$Sample_ID -hold_jid "hat01_"$Sample_ID "hatchet_pipe_02.count-alleles.sh" \
    # --TUMORS ${CASE_BAM_PATH} --SAMPLES ${ALLNAMES} --NORMAL ${CONTROL_BAM_PATH} --SNPS ${SNPS} --REF ${REF_hg38} \
    # --OUTPUTNORMAL ${OUTPUTNORMAL} --OUTPUTTUMORS ${OUTPUTTUMORS} --OUTPUTSNPS ${OUTPUTSNPS} \
    # --SAMTOOLS ${SAMTOOLS} --BCFTOOLS ${BCFTOOLS} --MINCOV ${MINCOV} --MAXCOV ${MAXCOV} --PROCESSES ${PROCESSES} --READQUALITY ${READQUALITY} --BASEQUALITY ${BASEQUALITY}



    #03. count-reads (WES에서 약 30분 걸림, /home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/hatchet/utils/)
    REFVERSION="hg38"  
    BAFFILE=${HATCHET_PATH_02}"/"${Sample_ID}"/"${Sample_ID}".outputnormal.tsv"
    OUTDIR=${HATCHET_PATH_03}"/"${Sample_ID}
    if [ -d ${OUTDIR} ]; then
        rm -rf ${OUTDIR}
    fi
    if [ ! -d ${OUTDIR} ]; then
        mkdir -p ${OUTDIR}
    fi
    for sublog in "hatchet_03.count-reads"; do
        if [ $logPath"/"$sublog ] ; then
            rm -rf $logPath"/"$sublog
        fi
        if [ ! -d $logPath"/"$sublog ] ; then
            mkdir -p $logPath"/"$sublog
        fi
    done
    qsub -pe smp 3 -o $logPath"/hatchet_03.count-reads" -e $logPath"/hatchet_03.count-reads" -N "hat03_"$Sample_ID -hold_jid "hat02_"$Sample_ID "hatchet_pipe_03.count-reads.sh" \
    --TUMORS ${CASE_BAM_PATH} --SAMPLES {$ALLNAMES} --NORMAL ${CONTROL_BAM_PATH} --BAFFILE ${BAFFILE} --REFVERSION ${REFVERSION} \
    --OUTDIR ${OUTDIR}  \
    --SAMTOOLS ${SAMTOOLS} --PROCESSES ${PROCESSES} --READQUALITY ${READQUALITY}



    #04. phase-snps  (이상하다)   python3 -m hatchet  download-panel -D "/home/goldpm1/Meningioma/11.cnv/3.hatchet/99.1000g" -R "1000GP_Phase3" )
    REFPANELDIR="/home/goldpm1/Meningioma/11.cnv/3.hatchet/99.1000g/"
    CHRNOTATION="True"
    OUTDIR=${HATCHET_PATH_04}"/"${Sample_ID}
    if [ ! -d ${OUTDIR} ]; then
        mkdir -p ${OUTDIR}
    fi
    # qsub -pe smp 3 -o $logPath"/hatchet_04.phase-snps" -e $logPath"/hatchet_04.phase-snps" -N "hat04_"$Sample_ID -hold_jid "hat03_"$Sample_ID "hatchet_pipe_04.phase-snps.sh" \
    # --SNPS ${SNPS} --REFPANELDIR ${REFPANELDIR} --REFGENOME ${REF_hg38} --REFVERSION ${REFVERSION} --CHRNOTATION ${CHRNOTATION} \
    # --OUTDIR ${OUTDIR}  \
    # --BCFTOOLS ${BCFTOOLS} --PROCESSES ${PROCESSES} --SHAPEIT ${SHAPEIT} --PICARD ${PICARD} --BGZIP ${BGZIP}


    #05. combine-counts
    ARRAY=${HATCHET_PATH_03}"/"${Sample_ID}
    BAFFILE=${OUTPUTTUMORS}
    TOTALCOUNTS=${HATCHET_PATH_03}"/"${Sample_ID}"/total.tsv"
    OUTFILE=${HATCHET_PATH_05}"/"${Sample_ID}"/"${Sample_ID}".tsv"
    MSR="5000"
    MTR="5000"
    for sublog in "hatchet_05.combine-counts"; do
        if [ $logPath"/"$sublog ] ; then
            rm -rf $logPath"/"$sublog
        fi
        if [ ! -d $logPath"/"$sublog ] ; then
            mkdir -p $logPath"/"$sublog
        fi
    done
    qsub -pe smp 3 -o $logPath"/hatchet_05.combine-counts" -e $logPath"/hatchet_05.combine-counts" -N "hat05_"$Sample_ID -hold_jid "hat03_"$Sample_ID",hat04_"$Sample_ID "hatchet_pipe_05.combine-counts.sh" \
    --ARRAY ${ARRAY} --BAFFILE ${BAFFILE} --TOTALCOUNTS ${TOTALCOUNTS} --REFVERSION ${REFVERSION} \
    --OUTFILE ${OUTFILE}  \
    --MSR ${MSR} --MTR ${MTR} --PROCESSES ${PROCESSES} 


    #06. cluster-bins
    INPUT_TSV=${OUTFILE}
    DECODING="map"
    OUTBINS=${HATCHET_PATH_06}"/"${Sample_ID}"/"${Sample_ID}".bbc"
    OUTSEGMENTS=${HATCHET_PATH_06}"/"${Sample_ID}"/"${Sample_ID}".seg"
    if [ ! -d ${OUTBINS%/*} ]; then
        mkdir -p ${OUTBINS%/*}
    fi
    for sublog in "hatchet_06.cluster-bins"; do
        if [ $logPath"/"$sublog ] ; then
            rm -rf $logPath"/"$sublog
        fi
        if [ ! -d $logPath"/"$sublog ] ; then
            mkdir -p $logPath"/"$sublog
        fi
    done
    qsub -pe smp 3 -o $logPath"/hatchet_06.cluster-bins" -e $logPath"/hatchet_06.cluster-bins" -N "hat06_"$Sample_ID -hold_jid "hat05_"$Sample_ID "hatchet_pipe_06.cluster-bins.sh" \
    --INPUT_TSV ${INPUT_TSV} --DECODING ${DECODING} --OUTBINS ${OUTBINS} --OUTSEGMENTS ${OUTSEGMENTS}

    #07. plot-bins
    INPUT_TSV=${OUTBINS}
    RUNDIR=${HATCHET_PATH_07}"/"${Sample_ID}
    if [ ! -d ${RUNDIR} ]; then
        mkdir -p ${RUNDIR}
    fi
    for sublog in "hatchet_07.plot-bins"; do
        if [ $logPath"/"$sublog ] ; then
            rm -rf $logPath"/"$sublog
        fi
        if [ ! -d $logPath"/"$sublog ] ; then
            mkdir -p $logPath"/"$sublog
        fi
    done
    qsub -pe smp 1 -o $logPath"/hatchet_07.plot-bins" -e $logPath"/hatchet_07.plot-bins" -N "hat07_"$Sample_ID -hold_jid "hat06_"$Sample_ID "hatchet_pipe_07.plot-bins.sh" \
    --INPUT_TSV ${INPUT_TSV} --RUNDIR ${RUNDIR}

    #08. compute_cnv   (/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/hatchet/bin/HATCHet.py)
    INPUT_BINS_SEGMENTS=${HATCHET_PATH_06}"/"${Sample_ID}"/"${Sample_ID}
    RUNNINGDIR=${HATCHET_PATH_08}"/"${Sample_ID}
    if [ ! -d ${RUNNINGDIR} ]; then
        mkdir -p ${RUNNINGDIR}
    fi
    for sublog in "hatchet_08.compute-cn"; do
        if [ $logPath"/"$sublog ] ; then
            rm -rf $logPath"/"$sublog
        fi
        if [ ! -d $logPath"/"$sublog ] ; then
            mkdir -p $logPath"/"$sublog
        fi
    done
    # qsub -pe smp 3 -o $logPath"/hatchet_08.compute-cn" -e $logPath"/hatchet_08.compute-cn" -N "hat08_"$Sample_ID -hold_jid "hat07_"$Sample_ID "hatchet_pipe_08.compute-cn.sh" \
    # --INPUT_BINS_SEGMENTS ${INPUT_BINS_SEGMENTS} --RUNNINGDIR ${RUNNINGDIR}
    

    #09. plot-cn

    #10. check


done


