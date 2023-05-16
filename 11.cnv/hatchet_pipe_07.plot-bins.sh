#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#export GRB_LICENSE_FILE=/opt/Yonsei/Gurobi/9.1.1/linux64/gurobi_nodes.lic

if ! options=$(getopt -o h --long INPUT_TSV:,RUNDIR:, -- "$@")
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
        --RUNDIR)
            RUNDIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


date

echo -e "python3 -m hatchet plot-bins  \
    ${INPUT_TSV} --rundir ${RUNDIR}  -c BAF -m Set1 --figsize 9,2.5 --markersize 1 --ymax 1 --ymin 0"


python3 -m hatchet plot-bins ${INPUT_TSV}  --rundir ${RUNDIR}  -c BAF -m Set1 --figsize 9,2.5 --markersize 50 --fontscale 2 --ymax 1 --ymin 0
python3 -m hatchet plot-bins ${INPUT_TSV}  --rundir ${RUNDIR} -c CBAF -m Set1 --figsize 9,2.5 --markersize 50 --fontscale 2 --ymax 1 --ymin 0

python3 -m hatchet plot-bins ${INPUT_TSV}  --rundir ${RUNDIR} -c BB --figsize 8,6 -m Set1 --markersize 50 --fontscale 2
python3 -m hatchet plot-bins ${INPUT_TSV}  --rundir ${RUNDIR} -c CBB  -m Set1 --markersize 50 --fontscale 2  --colwrap 1  -tS 0.005

python3 -m hatchet plot-bins ${INPUT_TSV}  --rundir ${RUNDIR} -c CRD -m Set1 --figsize 9,2.5 --markersize 30 --fontscale 2  # Read Depth ratio
