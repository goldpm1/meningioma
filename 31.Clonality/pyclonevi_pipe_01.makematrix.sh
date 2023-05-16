#!/bin/bash
#$ -S /bin/bash
#$ -cwd
# Basic Argument

if ! options=$(getopt -o h --long Sample_ID:,TISSUE:,SEQUENZA_MUTATION_PATH:,SEQUENZA_SEGMENT_PATH:,SEQUENZA_PLOIDY_PATH:,SEQUENZA_PURITY_PLOIDY_PATH:,SEQUENZA_TO_PYCLONEVI_MATRIX_PATH:,FACETCNV_TO_BED_DF_PATH:,FACETCNV_OUTPUT_PATH:,FACETCNV_PURITY_PLODY_PATH:,FACETCNV_TO_PYCLONEVI_MATRIX_PATH:,MUTECT_OUTPUT_PATH:,HC_OUTPUT_PATH:,HC_BLOOD_RANDOM_PICK_PATH:,GVCF_OUTPUT_PATH:, -- "$@")
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
        --Sample_ID)
            Sample_ID=$2
        shift 2 ;;
        --TISSUE)
            TISSUE=$2
        shift 2 ;;
        --SEQUENZA_MUTATION_PATH)
            SEQUENZA_MUTATION_PATH=$2
        shift 2 ;;
        --SEQUENZA_SEGMENT_PATH)
            SEQUENZA_SEGMENT_PATH=$2
        shift 2 ;;
        --SEQUENZA_PLOIDY_PATH)
            SEQUENZA_PLOIDY_PATH=$2
        shift 2 ;;
        --SEQUENZA_PURITY_PLOIDY_PATH)
            SEQUENZA_PURITY_PLOIDY_PATH=$2
        shift 2 ;;
        --SEQUENZA_TO_PYCLONEVI_MATRIX_PATH)
            SEQUENZA_TO_PYCLONEVI_MATRIX_PATH=$2
        shift 2 ;;
        --FACETCNV_OUTPUT_PATH)
            FACETCNV_OUTPUT_PATH=$2
        shift 2 ;;        
        --FACETCNV_TO_BED_DF_PATH)
            FACETCNV_TO_BED_DF_PATH=$2
        shift 2 ;;
        --FACETCNV_PURITY_PLODY_PATH)
            FACETCNV_PURITY_PLODY_PATH=$2
        shift 2 ;;
        --FACETCNV_TO_PYCLONEVI_MATRIX_PATH)
            FACETCNV_TO_PYCLONEVI_MATRIX_PATH=$2
        shift 2 ;;
        --MUTECT_OUTPUT_PATH)
            MUTECT_OUTPUT_PATH=$2
        shift 2 ;;
        --HC_OUTPUT_PATH)
            HC_OUTPUT_PATH=$2
        shift 2 ;;
        --HC_BLOOD_RANDOM_PICK_PATH)
            HC_BLOOD_RANDOM_PICK_PATH=$2
        shift 2 ;;
        --GVCF_OUTPUT_PATH)
            GVCF_OUTPUT_PATH=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



python3 pyclonevi_pipe_01.makematrix_sequenza_to_pyclonevi.py \
    --SEQUENZA_MUTATION_PATH ${SEQUENZA_MUTATION_PATH} \
    --SEQUENZA_SEGMENT_PATH ${SEQUENZA_SEGMENT_PATH}  \
    --SEQUENZA_PLOIDY_PATH ${SEQUENZA_PLOIDY_PATH} \
    --SEQUENZA_PURITY_PLOIDY_PATH ${SEQUENZA_PURITY_PLOIDY_PATH} \
    --SEQUENZA_TO_PYCLONEVI_MATRIX_PATH ${SEQUENZA_TO_PYCLONEVI_MATRIX_PATH} \
    --MUTECT_OUTPUT_PATH ${MUTECT_OUTPUT_PATH} \
    --HC_OUTPUT_PATH ${HC_OUTPUT_PATH} \
    --HC_BLOOD_RANDOM_PICK_PATH ${HC_BLOOD_RANDOM_PICK_PATH} \
    --GVCF_OUTPUT_PATH ${GVCF_OUTPUT_PATH} \
    --TISSUE ${TISSUE} \
    --Sample_ID ${Sample_ID}

python3 pyclonevi_pipe_01.makematrix_facetcnv_to_pyclonevi.py \
    --FACETCNV_OUTPUT_PATH ${FACETCNV_OUTPUT_PATH} \
    --FACETCNV_TO_BED_DF_PATH ${FACETCNV_TO_BED_DF_PATH} \
    --FACETCNV_TO_PYCLONEVI_MATRIX_PATH ${FACETCNV_TO_PYCLONEVI_MATRIX_PATH} \
    --FACETCNV_PURITY_PLODY_PATH ${FACETCNV_PURITY_PLODY_PATH} \
    --MUTECT_OUTPUT_PATH ${MUTECT_OUTPUT_PATH} \
    --HC_OUTPUT_PATH ${HC_OUTPUT_PATH} \
    --HC_BLOOD_RANDOM_PICK_PATH ${HC_BLOOD_RANDOM_PICK_PATH} \
    --GVCF_OUTPUT_PATH ${GVCF_OUTPUT_PATH}
    --TISSUE ${TISSUE} \
    --Sample_ID ${Sample_ID}