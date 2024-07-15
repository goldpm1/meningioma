#!/bin/bash
#$ -cwd
#$ -S /bin/bash

DIR=$1
SAMPLE=$2
REF=$3
idx=/data/project/linked_read/raw_data/genome.fa.idx
AnalysisPath=$4
mpileupPath=$5
BQ=$6
PROJECT=$7
SCRIPT=/data/project/linked_read/script

if ! options=$(getopt -o h --long REF:,BQ:,OUTPUT_MPILEUP_CALL_PATH:,OUTPUT_MPILEUP_PATH:,SAMPLE:, -- "$@")
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
        --REF)
            REF=$2
        shift 2 ;;
        --BQ)
            BQ=$2
        shift 2 ;;
        --OUTPUT_MPILEUP_CALL_PATH)
            OUTPUT_MPILEUP_CALL_PATH=$2
        shift 2 ;;
        --OUTPUT_MPILEUP_PATH)
            OUTPUT_MPILEUP_PATH=$2
        shift 2 ;;
        --SAMPLE)
            SAMPLE=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


# print start time
date

# 1. Checking Mutation for mpileup
cd $AnalysisPath
python3  "04.mpileup_01.checkmut_mq_indel.py" -q ${BQ} -F 1 -f 0 -o ${OUTPUT_MPILEUP_CALL_PATH} ${OUTPUT_MPILEUP_PATH}

# 2. Annotation
python3 "04.mpileup_02.parseCall.py" $AnalysisPath ${OUTPUT_MPILEUP_CALL_PATH} $SAMPLE ${REF} ${idx}

#python $SCRIPT/pysams/parsing_mpileup.py $DIR $SAMPLE 'B'$BQ $PROJECT

# print end time
date