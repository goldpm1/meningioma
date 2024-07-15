
#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long INPUT_BAM:,PRIMER_XML:,PRIMER:FASTA:,OUTPUT_BAM:, -- "$@")
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
        --PRIMER_XML)
            PRIMER_XML=$2
        shift 2 ;;
        --PRIMER_FASTA)
            PRIMER_FASTA=$2
        shift 2 ;;
        --INPUT_BAM)
            INPUT_BAM=$2
        shift 2 ;;
        --OUTPUT_BAM)
            OUTPUT_BAM=$2
        shift 2 ;;
        --)
            shift
            break
    esac
donec


echo -e "/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/isoseq refine ${INPUT_BAM} ${PRIMER_FASTA} ${OUTPUT_BAM}  --require-polya"

/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/isoseq refine ${INPUT_BAM} ${PRIMER_FASTA} ${OUTPUT_BAM} --require-polya