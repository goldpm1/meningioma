#!/bin/bash
#$ -S /bin/bash
#$ -cwd


if ! options=$(getopt -o h --long CASE_BAM:,CONTROL_BAM:,REF:,WIGGZ:,SEQUENZA_SEQZ:,SEQUENZA_SMALL_SEQZ:, -- "$@")
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
        --CASE_BAM)
            CASE_BAM=$2
        shift 2 ;;
        --CONTROL_BAM)
            CONTROL_BAM=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --WIGGZ)
            WIGGZ=$2
        shift 2 ;;
        --SEQUENZA_SEQZ)
            SEQUENZA_SEQZ=$2
        shift 2 ;;
        --SEQUENZA_SMALL_SEQZ)
            SEQUENZA_SMALL_SEQZ=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


## Print Start time
date

#1. Process a FASTA file to produce a GC Wiggle track file  (지금은 할 필요 없다)
# sequenza-utils gc_wiggle \
# -w 50 \
# --fasta $REF \
# -o ${WIGGZ}

#2. Process a BAM and Wiggle files to produce a seqz file
# sequenza-utils bam2seqz \
# -t ${CASE_BAM}  \
# -n ${CONTROL_BAM} \
# --fasta ${REF} \
# -gc ${WIGGZ} \
# -o ${SEQUENZA_SEQZ}

#3. Post-Process by binning the original seqz file
sequenza-utils seqz_binning \
--seqz ${SEQUENZA_SEQZ} \
-o ${SEQUENZA_SMALL_SEQZ} \
-w 50 \
--tabix /opt/Yonsei/htslib/1.14/bin/tabix


#4. chrM 이후로는 날려버리기

echo -e "python3 "sequenza_pipe_02.removechrM.py" \
    --INPUT_VCF_GZ ${SEQUENZA_SMALL_SEQZ} \
    --OUTPUT_VCF ${SEQUENZA_SMALL_SEQZ%.gz}"

python3 "sequenza_pipe_02.removechrM.py" \
    --INPUT_VCF_GZ ${SEQUENZA_SMALL_SEQZ} \
    --OUTPUT_VCF ${SEQUENZA_SMALL_SEQZ%.gz}


# Sort 하기  (맨 첫줄은 header니까 제외하고)
head -n +1 ${SEQUENZA_SMALL_SEQZ%.gz} > ${SEQUENZA_SMALL_SEQZ%.gz}".header" 
tail -n +2 ${SEQUENZA_SMALL_SEQZ%.gz} > ${SEQUENZA_SMALL_SEQZ%.gz}".unsorted.body"
sort -V -k1,1 -k2,2 ${SEQUENZA_SMALL_SEQZ%.gz}".unsorted.body" > ${SEQUENZA_SMALL_SEQZ%.gz}".sorted.body"
cat ${SEQUENZA_SMALL_SEQZ%.gz}".header" ${SEQUENZA_SMALL_SEQZ%.gz}".sorted.body" > ${SEQUENZA_SMALL_SEQZ%.gz}".sorted"


bgzip -c -f ${SEQUENZA_SMALL_SEQZ%.gz}".sorted" > ${SEQUENZA_SMALL_SEQZ}".removechrM.gz"
tabix  -f -s 1 -b 2 -e 2 -S 1 ${SEQUENZA_SMALL_SEQZ}".removechrM.gz"

rm -rf ${SEQUENZA_SMALL_SEQZ%.gz}".header"  ${SEQUENZA_SMALL_SEQZ%.gz}".unsorted.body" ${SEQUENZA_SMALL_SEQZ%.gz}".sorted.body" ${SEQUENZA_SMALL_SEQZ%.gz}".sorted"

date
