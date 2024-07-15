#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long REF:,INPUT_VCF:,INPUT_SNP_VCF:,INPUT_INDEL_VCF:,OUTPUT_SNP_VCF:,OUTPUT_INDEL_VCF:,OUTPUT_VCF:,TRANCHES_SNP:,RECAL_SNP:,TRANCHES_INDEL:,RECAL_INDEL:, -- "$@")
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
        --REF)
            REF=$2
        shift 2 ;;
        --INPUT_VCF)
            INPUT_VCF=$2
        shift 2 ;;
        --INPUT_SNP_VCF)
            INPUT_SNP_VCF=$2
        shift 2 ;;
        --INPUT_INDEL_VCF)
            INPUT_INDEL_VCF=$2
        shift 2 ;;
        --OUTPUT_SNP_VCF)
            OUTPUT_SNP_VCF=$2
        shift 2 ;;
        --OUTPUT_INDEL_VCF)
            OUTPUT_INDEL_VCF=$2
        shift 2 ;;
        --OUTPUT_VCF)
            OUTPUT_VCF=$2
        shift 2 ;;
        --TRANCHES_SNP)
            TRANCHES_SNP=$2
        shift 2 ;;
        --RECAL_SNP)
            RECAL_SNP=$2
        shift 2 ;;
        --TRANCHES_INDEL)
            TRANCHES_INDEL=$2
        shift 2 ;;
        --RECAL_INDEL)
            RECAL_INDEL=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


# gatk --java-options "-Xmx48g" SelectVariants \
#     -R ${REF} \
#     -V ${INPUT_VCF} \
#     --select-type-to-include SNP \
#     -O ${INPUT_SNP_VCF}

# gatk --java-options "-Xmx48g" SelectVariants \
#     -R ${REF} \
#     -V ${INPUT_VCF} \
#     --select-type-to-include INDEL \
#     -O ${INPUT_INDEL_VCF}


#### [SNP] ######


echo -e "gatk VariantRecalibrator \
   -R ${REF} \
   -V ${INPUT_SNP_VCF} \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /home/goldpm1/GATKresources/hapmap_3.3.hg38.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 /home/goldpm1/GATKresources/1000G_omni2.5.hg38.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /home/goldpm1/GATKresources/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /home/goldpm1/GATKresources/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   --max-gaussians 4 \
   --tranches-file ${TRANCHES_SNP} \
   -O ${RECAL_SNP}"

gatk VariantRecalibrator \
   -R ${REF} \
   -V ${INPUT_SNP_VCF} \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /home/goldpm1/GATKresources/hapmap_3.3.hg38.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 /home/goldpm1/GATKresources/1000G_omni2.5.hg38.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /home/goldpm1/GATKresources/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /home/goldpm1/GATKresources/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   --max-gaussians 2 \
   --tranches-file ${TRANCHES_SNP} \
   -O ${RECAL_SNP}
echo -e "VariantRecalibrator for SNP done"

 gatk ApplyVQSR \
     -R ${REF} \
     -V ${INPUT_SNP_VCF} \
     -O ${OUTPUT_SNP_VCF} \
     --truth-sensitivity-filter-level 98.0 \
     --tranches-file ${TRANCHES_SNP} \
     --recal-file ${RECAL_SNP}  \
     -mode SNP \
     --exclude-filtered true
echo -e "ApplyVQSR for SNP done"


# sort lexicographically
grep "^#" ${OUTPUT_SNP_VCF} > ${OUTPUT_SNP_VCF%vcf}"sorted.vcf" && grep -v "^#" ${OUTPUT_SNP_VCF} | sort -k1,1V -k2n >> ${OUTPUT_SNP_VCF%vcf}"sorted.vcf"
mv -f ${OUTPUT_SNP_VCF%vcf}"sorted.vcf" ${OUTPUT_SNP_VCF}



# ##### [INDEL] ######
# gatk VariantRecalibrator \
#     -R ${REF} \
#     -V ${INPUT_INDEL_VCF} \
#     --resource:mills,known=false,training=true,truth=true,prior=12.0 /home/goldpm1/GATKresources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
#     --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /home/goldpm1/GATKresources/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
#     -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
#     --mode INDEL \
#     --max-gaussians 2 \
#     --tranches-file ${TRANCHES_INDEL} \
#     -O ${RECAL_INDEL}
#     #-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0  \

#  gatk ApplyVQSR \
#      -R ${REF} \
#      -V ${INPUT_INDEL_VCF} \
#      -O ${OUTPUT_INDEL_VCF} \
#      --truth-sensitivity-filter-level 90.0 \
#      --tranches-file ${TRANCHES_INDEL} \
#      --recal-file ${RECAL_INDEL}  \
#      -mode INDEL \
#      --exclude-filtered true




rm -rf ${TRANCHES_SNP} ${TRANCHES_SNP}".idx" ${RECAL_SNP} ${RECAL_SNP}".idx" ${TRANCHES_INDEL} ${TRANCHES_INDEL}".idx" ${RECAL_INDEL} ${RECAL_INDEL}".idx"



# 하나로 합쳐주기
gatk  MergeVcfs -I ${OUTPUT_SNP_VCF} -I ${INPUT_INDEL_VCF} -O ${OUTPUT_VCF}
bgzip -c -f ${OUTPUT_VCF} > ${OUTPUT_VCF}".gz"
tabix -p vcf ${OUTPUT_VCF}".gz"




# # 왜인지 모르겠으나 2개 중복이 생겨서 하나를 제거하는 과정

# # awk '!seen[$0]++' $Output"/"${Filename%.vcf*}".VQSR.vcf" > $Output"/"${Filename%.vcf*}".VQSR2.vcf"
# # grep '#' $Output"/"${Filename%.vcf*}".VQSR2.vcf" > $Output"/"${Filename%.vcf*}".VQSR.vcf" && \
# # grep "PASS" $Output"/"${Filename%.vcf*}".VQSR2.vcf" | grep "GT:AD:DP" | grep -v "Description" >> $Output"/"${Filename%.vcf*}".VQSR.vcf"
# # rm -rf $Output"/"${Filename%.vcf*}".VQSR2.vcf"

# # bgzip -c $Output"/"${Filename%.vcf*}".VQSR.vcf" > $Output"/"${Filename%.vcf*}".VQSR.vcf.gz"
# # tabix -p vcf $Output"/"${Filename%.vcf*}".VQSR.vcf.gz"

# # rm -rf $Output"/"${Filename%.vcf*}".VQSR2.vcf"