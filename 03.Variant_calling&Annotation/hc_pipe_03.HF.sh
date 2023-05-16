#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if ! options=$(getopt -o h --long INPUT_VCF_GZ:,OUTPUT_VCF_GZ:,OUTPUT_VCF:,REF:, -- "$@")
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
        --INPUT_VCF_GZ)
            INPUT_VCF_GZ=$2
        shift 2 ;;
        --OUTPUT_VCF_GZ)
            OUTPUT_VCF_GZ=$2
        shift 2 ;;
        --OUTPUT_VCF)
            OUTPUT_VCF=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

rm -rf ${OUTPUT_VCF} ${OUTPUT_VCF_GZ} ${OUTPUT_VCF_GZ}".tbi" ${OUTPUT_VCF%vcf}"temp.vcf.gz" ${OUTPUT_VCF%vcf}"temp.vcf.idx"

# 그냥 hard filter하는게 마음 편하다.  DP > 100, Alt > 5, Ref > 5

python3 hc_pipe_hardfilter.py \
    --INPUT_VCF ${INPUT_VCF_GZ} \
    --OUTPUT_VCF ${OUTPUT_VCF} \
    --SAMPLE_THRESHOLD "Blood,Tumor,Dura" \
    --DP_THRESHOLD 100 --REF_THRESHOLD 5 --ALT_THRESHOLD 5 --REMOVE_MULTIALLELIC True --PASS False --REMOVE_MITOCHONDRIAL_DNA True \

bgzip -c -f ${OUTPUT_VCF} > ${OUTPUT_VCF_GZ} 
tabix -p vcf ${OUTPUT_VCF_GZ}



# gatk VariantFiltration \
#     -R ${REF} \
#     -V ${INPUT_VCF_GZ} \
#     -O ${OUTPUT_VCF%vcf}"temp.vcf" \
#     --filter-name "DPfilter" --filter-expression "DP>100" \
    # --filter-name "HOMOfilter" --filter-expression "AF=1.00"
# grep '#' ${OUTPUT_VCF%vcf}"temp.vcf" > ${OUTPUT_VCF}
# grep -v '#' ${OUTPUT_VCF%vcf}"temp.vcf" | grep 'PASS' >> ${OUTPUT_VCF}

# bgzip -c -f ${OUTPUT_VCF} > ${OUTPUT_VCF_GZ} 
# tabix -p vcf ${OUTPUT_VCF_GZ}

# rm -rf ${OUTPUT_VCF%vcf}"temp.vcf"
