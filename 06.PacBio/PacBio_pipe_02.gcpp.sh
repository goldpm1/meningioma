#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long BAM_PATH:,REF:,REF_MMI:,PON:,gnomad:,TMP_PATH:,OUTPUT_VCF_PATH:,OUTPUT_MUTECT_GZ:,OUTPUT_MUTECT:,OUTPUT_FMC_PATH:,OUTPUT_FMC_HF_PATH:,OUTPUT_FMC_HF_RMBLACK_PATH:,SAMPLE_THRESHOLD:,DP_THRESHOLD:,ALT_THRESHOLD:,REMOVE_MULTIALLELIC:,PASS:,REMOVE_MITOCHONDRIAL_DNA:,BLACKLIST:, -- "$@")
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
        --BAM_PATH)
            BAM_PATH=$2
        shift 2 ;;
        --REF)
            REF=$2
        shift 2 ;;
        --REF_MMI)
            REF_MMI=$2
        shift 2 ;;
        --PON)
            PON=$2
        shift 2 ;;
        --gnomad)
            gnomad=$2
        shift 2 ;;
        --TMP_PATH)
            TMP_PATH=$2
        shift 2 ;;
        --OUTPUT_VCF_PATH)
            OUTPUT_VCF_PATH=$2
        shift 2 ;;
        --OUTPUT_MUTECT_GZ)
            OUTPUT_MUTECT_GZ=$2
        shift 2 ;;
        --OUTPUT_MUTECT)
            OUTPUT_MUTECT=$2
        shift 2 ;;
        --OUTPUT_FMC_PATH)
            OUTPUT_FMC_PATH=$2
        shift 2 ;;
        --OUTPUT_FMC_HF_PATH)
            OUTPUT_FMC_HF_PATH=$2
        shift 2 ;;
        --OUTPUT_FMC_HF_RMBLACK_PATH)
            OUTPUT_FMC_HF_RMBLACK_PATH=$2
        shift 2 ;;
        --SAMPLE_THRESHOLD)
            SAMPLE_THRESHOLD=$2
        shift 2 ;;
        --DP_THRESHOLD)
            DP_THRESHOLD=$2
        shift 2 ;;
        --ALT_THRESHOLD)
            ALT_THRESHOLD=$2
        shift 2 ;;
        --REMOVE_MULTIALLELIC)
            REMOVE_MULTIALLELIC=$2
        shift 2 ;;
        --PASS)
            PASS=$2
        shift 2 ;;
        --REMOVE_MITOCHONDRIAL_DNA)
            REMOVE_MITOCHONDRIAL_DNA=$2
        shift 2 ;;
        --BLACKLIST)
            BLACKLIST=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

#1. gcpp
# /opt/Yonsei/pacbio/smrtlink/smrtcmds/bin/gcpp \
#     -r "/home/goldpm1/miniconda3/envs/PacBio/lib/python3.7/site-packages/pbcore/data/lambdaNEB.fa" \
#     -o ${OUTPUT_VCF_PATH}
#     ${BAM_PATH} \




# 2. Mutect
#1.Paired sample
gatk --java-options "-Xmx24g" Mutect2 \
-R ${REF} \
-I ${BAM_PATH} \
--panel-of-normals ${PON} \
--germline-resource ${gnomad} \
-O ${OUTPUT_MUTECT_GZ} \
--tmp-dir ${TMP_PATH} \
--dont-use-soft-clipped-bases true \
--intervals chr22:29000000-31000000


gunzip -f ${OUTPUT_MUTECT_GZ}

# 추후 MutationTimeR를 위해서 FORMAT 을 바꿔주기
sed 's/ID=AD,Number=R/ID=AD,Number=2/' ${OUTPUT_MUTECT} > ${OUTPUT_MUTECT}".temp"
mv ${OUTPUT_MUTECT}".temp" ${OUTPUT_MUTECT}
rm -rf ${OUTPUT_MUTECT}".temp" ${OUTPUT_MUTECT_GZ}

bgzip -c -f ${OUTPUT_MUTECT} > ${OUTPUT_MUTECT_GZ}
tabix -f -p vcf ${OUTPUT_MUTECT_GZ}


gatk FilterMutectCalls -R ${REF} -V ${OUTPUT_MUTECT_GZ} -O ${OUTPUT_FMC_PATH} --max-events-in-region 1 --min-median-read-position 8 --min-median-base-quality 20  --min-reads-per-strand 1 --min-median-mapping-quality 20

#3. Hard filter
python3 "/data/project/Meningioma/script/03.Variant_calling&Annotation/mutect_pipe_hardfilter.py" --INPUT_VCF ${OUTPUT_FMC_PATH} --OUTPUT_VCF ${OUTPUT_FMC_HF_PATH} --SAMPLE_THRESHOLD ${SAMPLE_THRESHOLD} --DP_THRESHOLD ${DP_THRESHOLD} --ALT_THRESHOLD ${ALT_THRESHOLD}  --REMOVE_MULTIALLELIC ${REMOVE_MULTIALLELIC} --PASS ${PASS} --REMOVE_MITOCHONDRIAL_DNA ${REMOVE_MITOCHONDRIAL_DNA}

#4. Repeat region 지우기
bedtools intersect -header -v -a ${OUTPUT_FMC_HF_PATH} -b ${BLACKLIST} > ${OUTPUT_FMC_HF_RMBLACK_PATH}


# bgzip & tabix
bgzip -c -f ${OUTPUT_FMC_HF_RMBLACK_PATH} > ${OUTPUT_FMC_HF_RMBLACK_PATH}".gz"
tabix -p vcf ${OUTPUT_FMC_HF_RMBLACK_PATH}".gz"
cp ${OUTPUT_FMC_HF_RMBLACK_PATH}".gz.tbi" ${OUTPUT_FMC_HF_RMBLACK_PATH}".tbi"



# VEP
bash "/data/project/Meningioma/script/03.Variant_calling&Annotation/mutect_pair_pipe_20.vep.sh" \
    --REF ${REF} \
    --INPUT_VCF ${OUTPUT_FMC_HF_RMBLACK_PATH} \
    --OUTPUT_VCF ${OUTPUT_FMC_HF_RMBLACK_PATH%vcf}"vep.vcf"
