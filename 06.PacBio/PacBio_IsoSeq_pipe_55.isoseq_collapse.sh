
#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long INPUT_BAM:,OUTPUT_GFF:, -- "$@")
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
        --INPUT_BAM)
            INPUT_BAM=$2
        shift 2 ;;
        --OUTPUT_GFF)
            OUTPUT_GFF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



echo -e "/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/isoseq collapse ${INPUT_BAM} ${OUTPUT_GFF}"
/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/isoseq collapse ${INPUT_BAM} ${OUTPUT_GFF} 


# isoseq collapse --do-not-collapse-extra-5exons aligned.sorted.bam out.gff or
# isoseq collapse --do-not-collapse-extra-5exons aligned.sorted.bam ccs.bam out.gff