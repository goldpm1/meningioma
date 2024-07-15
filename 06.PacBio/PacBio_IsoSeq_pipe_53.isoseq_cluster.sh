
#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long FILE_LIST:,INPUT_BAM:,OUTPUT_BAM:, -- "$@")
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
        --FILE_LIST)
            FILE_LIST=$2
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
done



# echo -e "/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/isoseq cluster ${FILE_LIST} ${OUTPUT_BAM} --verbose --use-qvs"
# /opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/isoseq cluster ${FILE_LIST} ${OUTPUT_BAM} --verbose --use-qvs


echo -e "/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/isoseq cluster2 ${INPUT_BAM} ${OUTPUT_BAM}"
/opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/isoseq cluster2 ${INPUT_BAM} ${OUTPUT_BAM} 
