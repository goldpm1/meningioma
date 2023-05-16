#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#module load gurobi
#export GRB_LICENSE_FILE=/opt/Yonsei/Gurobi/9.1.1/linux64/gurobi_nodes.lic
export GRB_LICENSE_FILE=/opt/Yonsei/Gurobi/10.0.1/linux64/gurobi.lic


if ! options=$(getopt -o h --long CONTROL_BAM_PATH:,CASE_BAM_PATH:,SAMPLES:,ALLNAMES:,HATCHET_PATH:,REF:,REF_VERSION:,dbSNP:,REGION:,SAMTOOLS:,BCFTOOLS:,BGZIP:,SHAPEIT:,PICARD:,MINCOV:,MAXCOV:,MINREADS:,MAXREADS:,PROCESSES:,READQUALITY:,BASEQUALITY:,CHR_NOTATION:,BIN:,PHASE:,RANDOM:,LOGPATH:,SAMPLE_ID:, -- "$@")
then
    echo "ERROR: invalid options"
    exit 1
fi

eval set -- $options

while true; do
    case "$1" in
        -h|--help)
            echo "Usage"
        shift ;;
        --CONTROL_BAM_PATH)
            CONTROL_BAM_PATH=$2
        shift 2 ;;
        --CASE_BAM_PATH)
            CASE_BAM_PATH=$2
        shift 2 ;;
        --SAMPLES)
            SAMPLES=$2
        shift 2 ;;
        --ALLNAMES)
            ALLNAMES=$2
        shift 2 ;;
        --HATCHET_PATH)
            HATCHET_PATH=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --REF_VERSION)
            REF_VERSION=$2
        shift 2 ;;
        --dbSNP)
            dbSNP=$2
        shift 2 ;;
        --REGION)
            REGION=$2
        shift 2 ;;
        --SAMTOOLS)
            SAMTOOLS=$2
        shift 2 ;;
        --BCFTOOLS)
            BCFTOOLS=$2
        shift 2 ;;
        --BGZIP)
            BGZIP=$2
        shift 2 ;;
        --SHAPEIT)
            SHAPEIT=$2
        shift 2 ;;
        --PICARD)
            PICARD=$2
        shift 2 ;;
        --MINCOV)
            MINCOV=$2
        shift 2 ;;
        --MAXCOV)
            MAXCOV=$2
        shift 2 ;;
        --MINREADS)
            MINREADS=$2
        shift 2 ;;
        --MAXREADS)
            MAXREADS=$2
        shift 2 ;;
        --PROCESSES)
            PROCESSES=$2
        shift 2 ;;
        --READQUALITY)
            READQUALITY=$2
        shift 2 ;;
        --BASEQUALITY)
            BASEQUALITY=$2
        shift 2 ;;
        --CHR_NOTATION)
            CHR_NOTATION=$2
        shift 2 ;;
        --BIN)
            BIN=$2
        shift 2 ;;
        --PHASE)
            PHASE=$2
        shift 2 ;;
        --RANDOM)
            RANDOM=$2
        shift 2 ;;
        --LOGPATH)
            logPath=$2
        shift 2 ;;
        --SAMPLE_ID)
            Sample_ID=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



CASE_BAM_PATH_NEW=$(echo ${CASE_BAM_PATH} | sed 's/[{}]//g' | sed 's/,/ /g')
SAMPLES_NEW=$(echo ${SAMPLES} | sed 's/[{}]//g' | sed 's/,/ /g')
ALLNAMES_NEW=$(echo ${ALLNAMES} | sed 's/[{}]//g' | sed 's/,/ /g')

# make subpath
HATCHET_PATH_01=${HATCHET_PATH}"/01.RDR/"${Sample_ID}
HATCHET_PATH_02=${HATCHET_PATH}"/02.SNP/"${Sample_ID}
HATCHET_PATH_03=${HATCHET_PATH}"/03.BAF/"${Sample_ID}
HATCHET_PATH_04=${HATCHET_PATH}"/04.BB/"${Sample_ID}
HATCHET_PATH_05=${HATCHET_PATH}"/05.BBC/"${Sample_ID}
HATCHET_PATH_06=${HATCHET_PATH}"/06.PLO/"${Sample_ID}
HATCHET_PATH_07=${HATCHET_PATH}"/07.RES/"${Sample_ID}
HATCHET_PATH_08=${HATCHET_PATH}"/08.SUM/"${Sample_ID}

for subpath in ${HATCHET_PATH_01} ${HATCHET_PATH_02} ${HATCHET_PATH_03} ${HATCHET_PATH_04} ${HATCHET_PATH_05} ${HATCHET_PATH_06} ${HATCHET_PATH_07} ${HATCHET_PATH_08} ; do
    if [ ! -d ${subpath} ]; then
        mkdir -p ${subpath}
    fi
done

#01. binBAM (Normal, Tumor에서 BAM을 쪼갬)
date
# echo -e "\n/opt/Yonsei/python/3.8.1/bin/python3 -m hatchet binBAM -N ${CONTROL_BAM_PATH} -T ${CASE_BAM_PATH_NEW} -S ${ALLNAMES_NEW} -b ${BIN} -g ${REF} -r ${REGION} -j ${PROCESSES} -q ${READQUALITY} -O ${HATCHET_PATH_01}/normal.1bed -o ${HATCHET_PATH_01}/tumor.1bed -t ${HATCHET_PATH_01}/total.tsv"
# /opt/Yonsei/python/3.8.1/bin/python3 -m hatchet binBAM -N ${CONTROL_BAM_PATH} -T ${CASE_BAM_PATH_NEW} -S ${ALLNAMES_NEW} -b ${BIN} -g ${REF} -r ${REGION} -j ${PROCESSES} -q ${READQUALITY} -O ${HATCHET_PATH_01}"/normal.1bed" -o ${HATCHET_PATH_01}"/tumor.1bed" -t ${HATCHET_PATH_01}"/total.tsv"

#02. SNPcaller (Normal에서 germline hetero SNP를 call함)
# date
# echo -e "\n/opt/Yonsei/python/3.8.1/bin/python3 -m hatchet SNPCaller -N ${CONTROL_BAM_PATH} -r ${REF} -j ${PROCESSES} -c ${MINREADS} -C ${MAXREADS} -R ${dbSNP} -q ${READQUALITY} -Q ${BASEQUALITY} -o ${HATCHET_PATH_02}"
# /opt/Yonsei/python/3.8.1/bin/python3 -m hatchet SNPCaller -N ${CONTROL_BAM_PATH} -r ${REF} -j ${PROCESSES} -c ${MINREADS} -C ${MAXREADS} -R ${dbSNP} -q ${READQUALITY} -Q ${BASEQUALITY} -o ${HATCHET_PATH_02}

# #03. deBAF (count-alleles에 해당함.   여기도 마찬가지로 -l (OUTPUTSNPS)가 없다)
# date
# echo -e "\n/opt/Yonsei/python/3.8.1/bin/python3 -m hatchet deBAF -st ${SAMTOOLS} -bt ${BCFTOOLS} -N ${CONTROL_BAM_PATH} -T ${CASE_BAM_PATH_NEW} -S ${ALLNAMES_NEW} -r ${REF} -j ${PROCESSES} -L ${HATCHET_PATH_02}/*.vcf.gz -c ${MINREADS} -C ${MAXREADS} -O ${HATCHET_PATH_03}/normal.1bed -o ${HATCHET_PATH_03}/tumor.1bed -l "${HATCHET_PATH_03}
# /opt/Yonsei/python/3.8.1/bin/python3 -m hatchet deBAF -st ${SAMTOOLS} -bt ${BCFTOOLS} -N ${CONTROL_BAM_PATH} -T ${CASE_BAM_PATH_NEW} -S ${ALLNAMES_NEW} -r ${REF} -j ${PROCESSES} -L ${HATCHET_PATH_02}/*.vcf.gz -c ${MINREADS} -C ${MAXREADS} -O ${HATCHET_PATH_03}"/normal.1bed" -o ${HATCHET_PATH_03}"/tumor.1bed" -l ${HATCHET_PATH_03}

# #04. comBBO  (count-reads, combine-counts에 해당하는 것: BAF는 tumor BAF를 사용한다)
# date
# echo -e "\n/opt/Yonsei/python/3.8.1/bin/python3 -m hatchet comBBo -c ${HATCHET_PATH_01}/normal.1bed -C ${HATCHET_PATH_01}/tumor.1bed -B ${HATCHET_PATH_03}/tumor.1bed -t ${HATCHET_PATH_01}/total.tsv -p ${PHASE} > ${HATCHET_PATH_04}/bulk.bb"
# /opt/Yonsei/python/3.8.1/bin/python3 -m hatchet comBBo -c ${HATCHET_PATH_01}"/normal.1bed" -C ${HATCHET_PATH_01}"/tumor.1bed" -B ${HATCHET_PATH_03}"/tumor.1bed" -t ${HATCHET_PATH_01}"/total.tsv" -p ${PHASE} > ${HATCHET_PATH_04}"/bulk.bb"

# #05. cluBB
# date
# echo -e "\n/opt/Yonsei/python/3.8.1/bin/python3 -m hatchet cluBB ${HATCHET_PATH_04}/bulk.bb -o ${HATCHET_PATH_05}/bulk.seg -O ${HATCHET_PATH_05}/bulk.bbc -e ${RANDOM} -d 0.1 -tR 0.5 -tB 0.04"
# /opt/Yonsei/python/3.8.1/bin/python3 -m hatchet cluBB ${HATCHET_PATH_04}"/bulk.bb" -o ${HATCHET_PATH_05}"/bulk.seg" -O ${HATCHET_PATH_05}"/bulk.bbc" -e ${RANDOM} -d 0.1 -tR 0.5 -tB 0.04

# #06. Plot
# rm -rf ${HATCHET_PATH_06}
# mkdir ${HATCHET_PATH_06}
# date
# /opt/Yonsei/python/3.8.1/bin/python3 -m hatchet BBot  -c BAF ${HATCHET_PATH_05}"/bulk.bbc" --figsize 6,3 --markersize 20 --ymin 0 --ymax 1 --fontscale 2  --rundir ${HATCHET_PATH_06}
# /opt/Yonsei/python/3.8.1/bin/python3 -m hatchet BBot  -c CBAF ${HATCHET_PATH_05}"/bulk.bbc" --figsize 6,3 --markersize 20 --ymin 0 --ymax 1 --fontscale 2  --rundir ${HATCHET_PATH_06}
# /opt/Yonsei/python/3.8.1/bin/python3 -m hatchet BBot  -c RD   ${HATCHET_PATH_05}"/bulk.bbc" --figsize 6,3 --markersize 20  --ymin 0 --ymax 2  --fontscale 2 --rundir ${HATCHET_PATH_06}
# /opt/Yonsei/python/3.8.1/bin/python3 -m hatchet BBot  -c BB    ${HATCHET_PATH_05}"/bulk.bbc" --markersize 20  --fontscale 2 --rundir ${HATCHET_PATH_06}
# /opt/Yonsei/python/3.8.1/bin/python3 -m hatchet BBot  -c CBB ${HATCHET_PATH_05}"/bulk.bbc" --markersize 20  --fontscale 2 -tS 0.005 --ymin 0 --ymax 4 --colwrap 1 --rundir ${HATCHET_PATH_06}


# cd ../${RES} (compute-cn과 같은 기능)
echo -e "/opt/Yonsei/python/3.8.1/bin/python3 -m hatchet solve -i ${HATCHET_PATH_05}/bulk -n 2,6 -p 400 -u 0.06 -eD 6 -eT 12 -g 0.35 -l 0.6 -j ${PROCESSES} -r ${RANDOM} --diploid --runningdir ${HATCHET_PATH_07}"
/opt/Yonsei/python/3.8.1/bin/python3 -m hatchet solve -i ${HATCHET_PATH_05}"/bulk" -n 2,6 -p 400 -u 0.06 -eD 6 -eT 12 -g 0.35 -l 0.6 -j ${PROCESSES} -r ${RANDOM} --diploid --runningdir ${HATCHET_PATH_07}

# cd ../${SUM}
# /opt/Yonsei/python/3.8.1/bin/python3 -m hatchet BBeval ../${RES}/best.bbc.ucn -rC 10 -rG 1



