#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long TISSUE:,BAMTYPE:, -- "$@")
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
        --TISSUE)
            TISSUE=$2
        shift 2 ;;
        --BAMTYPE)
            BAMTYPE=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e  "python3 03.Amplicon_matrix_by_pysam.py --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}"
python3 "03.Amplicon_matrix_by_pysam.py" --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}
