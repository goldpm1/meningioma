#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#export GRB_LICENSE_FILE=/opt/Yonsei/Gurobi/9.1.1/linux64/gurobi_nodes.lic

if ! options=$(getopt -o h --long SNPS:,REFPANELDIR:,REFGENOME:,REFVERSION:,CHRNOTATION:,OUTDIR:,BCFTOOLS:,SHAPEIT:,PICARD:,BGZIP:,PROCESSES:, -- "$@")
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
        --SNPS)
            SNPS=$2
        shift 2 ;;
        --REFPANELDIR)
            REFPANELDIR=$2
        shift 2 ;;
        --REFGENOME)
            REFGENOME=$2
        shift 2 ;;
        --REFVERSION)
            REFVERSION=$2
        shift 2 ;;
        --CHRNOTATION)
            CHRNOTATION=$2
        shift 2 ;;
        --OUTDIR)
            OUTDIR=$2
        shift 2 ;;
        --BCFTOOLS)
            BCFTOOLS=$2
        shift 2 ;;
        --PROCESSES)
            PROCESSES=$2
        shift 2 ;;
        --SHAPEIT)
            SHAPEIT=$2
        shift 2 ;;
        --PICARD)
            PICARD=$2
        shift 2 ;;
        --BGZIP)
            BGZIP=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


date

SNPS_NEW=${SNPS}"/*.vcf.gz"

echo -e "python3 -m hatchet phase-snps \
    --snps ${SNPS_NEW} --refpaneldir ${REFPANELDIR} --refgenome ${REFGENOME} --refversion ${REFVERSION} -N ${CHRNOTATION} \
    --outdir ${OUTDIR} \
    --bcftools ${BCFTOOLS} --processes ${PROCESSES}  --shapeit ${SHAPEIT} --picard ${PICARD} --bgzip ${BGZIP}"

python3 -m hatchet phase-snps \
    --snps ${SNPS_NEW} --refpaneldir ${REFPANELDIR} --refgenome ${REFGENOME} --refversion ${REFVERSION} -N ${CHRNOTATION} \
    --outdir ${OUTDIR} \
    --bcftools ${BCFTOOLS} --processes ${PROCESSES}  --shapeit ${SHAPEIT} --picard ${PICARD} --bgzip ${BGZIP}
