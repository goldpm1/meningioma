#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if ! options=$(getopt -o h --long INPUT_VCF:,RECAL_FILE:,TRANCHES_FILE:,OUTPUT_VCF_GZ:,OUTPUT_VCF:,REF:, -- "$@")
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
        --INPUT_VCF)
            INPUT_VCF=$2
        shift 2 ;;
        --RECAL_FILE)
            RECAL_FILE=$2
        shift 2 ;;
        --TRANCHES_FILE)
            TRANCHES_FILE=$2
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


gatk VariantRecalibrator \
   -R ${REF} \
   -V ${INPUT_VCF} \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /home/goldpm1/GATKresources/hapmap_3.3.hg38.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 /home/goldpm1/GATKresources/1000G_omni2.5.hg38.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /home/goldpm1/GATKresources/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /home/goldpm1/GATKresources/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode BOTH \
   -O ${RECAL_FILE} \
   --tranches-file ${TRANCHES_FILE}


 gatk ApplyVQSR \
     -R ${REF} \
     -V ${INPUT_VCF} \
     -O ${OUTPUT_VCF_GZ} \
     --truth-sensitivity-filter-level 99.0 \
     --tranches-file ${TRANCHES_FILE} \
     --recal-file ${RECAL_FILE}  \
     -mode BOTH \
     --exclude-filtered true

  gunzip ${OUTPUT_VCF_GZ} 
  bgzip -c -f ${OUTPUT_VCF} > ${OUTPUT_VCF_GZ}
  tabix -p vcf ${OUTPUT_VCF_GZ}
