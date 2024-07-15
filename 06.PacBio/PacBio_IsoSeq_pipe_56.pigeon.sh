#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long GENCODE_GTF:,REF:,ID:,INPUT_GFF:,OUTPUT_DIR:, -- "$@")
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
        --GENCODE_GTF)
            GENCODE_GTF=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --ID)
            ID=$2
        shift 2 ;;
        --INPUT_GFF)
            INPUT_GFF=$2
        shift 2 ;;
        --OUTPUT_DIR)
            OUTPUT_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e "/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pigeon prepare ${GENCODE_GTF} ${REF}"
echo -e "/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pigeon prepare ${INPUT_GFF}"
echo -e "/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pigeon classify ${INPUT_GFF%gff}"sorted.gff" ${GENCODE_GTF%gtf}"sorted.gtf" ${REF}"


# /opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pigeon prepare ${GENCODE_GTF} ${REF} 
# /opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pigeon prepare ${INPUT_GFF}
#/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pigeon classify ${INPUT_GFF%gff}"sorted.gff" ${GENCODE_GTF%gtf}"sorted.gtf" ${REF} --fl ${INPUT_GFF%gff}"flnc_count.txt" -d ${OUTPUT_DIR}
/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/pigeon filter ${OUTPUT_DIR}"/"${ID}"_classification.txt" --isoforms ${INPUT_GFF%gff}"sorted.gff" --mono-exon
