#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long Sample_ID:,SNV_INPUT:,CNV_INPUT:,OUTPUT_PATH:,OUTPUT_PLOT_DIR:,MIN_VAF:,N_cutoff:,pi_cutoff:, -- "$@")
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
        --SNV_INPUT)
            SNV_INPUT=$2
        shift 2 ;;
        --CNV_INPUT)
            CNV_INPUT=$2
        shift 2 ;;
        --Sample_ID)
            Sample_ID=$2
        shift 2 ;;
        --OUTPUT_PATH)
            OUTPUT_PATH=$2
        shift 2 ;;
        --OUTPUT_PLOT_DIR)
            OUTPUT_PLOT_DIR=$2
        shift 2 ;;
        --MIN_VAF)
            MIN_VAF=$2
        shift 2 ;;
        --N_cutoff)
            N_cutoff=$2
        shift 2 ;;
        --pi_cutoff)
            pi_cutoff=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



echo "Start time: $(date)"
eval "$(mamba shell hook --shell bash)"
mamba activate TINC

echo "Rscript TINC_pipe_02.run.R --Sample_ID ${Sample_ID} --SNV_INPUT ${SNV_INPUT} --CNV_INPUT ${CNV_INPUT} --OUTPUT_PATH ${OUTPUT_PATH} --OUTPUT_PLOT_DIR ${OUTPUT_PLOT_DIR} --MIN_VAF ${MIN_VAF} --N_cutoff ${N_cutoff} --pi_cutoff ${pi_cutoff}"
Rscript TINC_pipe_02.run.R --Sample_ID ${Sample_ID} --SNV_INPUT ${SNV_INPUT} --CNV_INPUT ${CNV_INPUT} --OUTPUT_PATH ${OUTPUT_PATH} --OUTPUT_PLOT_DIR ${OUTPUT_PLOT_DIR} --MIN_VAF ${MIN_VAF} --N_cutoff ${N_cutoff} --pi_cutoff ${pi_cutoff}

mamba deactivate
echo "End time: $(date)"