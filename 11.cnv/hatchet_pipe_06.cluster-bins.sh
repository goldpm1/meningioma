#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#export GRB_LICENSE_FILE=/opt/Yonsei/Gurobi/9.1.1/linux64/gurobi_nodes.lic

if ! options=$(getopt -o h --long INPUT_TSV:,DECODING:,OUTBINS:,OUTSEGMENTS:, -- "$@")
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
        --INPUT_TSV)
            INPUT_TSV=$2
        shift 2 ;;
        --DECODING)
            DECODING=$2
        shift 2 ;;
        --OUTBINS)
            OUTBINS=$2
        shift 2 ;;
        --OUTSEGMENTS)
            OUTSEGMENTS=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


date


echo -e "python3 -m hatchet cluster-bins ${INPUT_TSV} \
    --decoding ${DECODING} --outbins ${OUTBINS} --outsegments ${OUTSEGMENTS}"

python3 -m hatchet cluster-bins ${INPUT_TSV} \
    --decoding ${DECODING} --outbins ${OUTBINS} --outsegments ${OUTSEGMENTS}
