#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long INPUT_VCF:,REF:,OUTPUT_VCF:, -- "$@")
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
        --OUTPUT_VCF)
            OUTPUT_VCF=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

Resource="/data/project/DC_WGS/Resource/"

####### 101.0을 쓰는 것을 추천
# /opt/Yonsei/ensembl-vep/104.3/vep \
# -v -assembly "GRCh38" \
# --terms "SO" \
# --fasta $REF \
# --no_stats \
# --everything \
# -i ${INPUT_VCF} --format "vcf" \
# -o ${OUTPUT_VCF} --force \
# --dir_plugins /data/public/VEP/104/Plugins \
# --custom ${Resource}dbSNP/GCF_000001405.39.re.vep_custom2.vcf.gz,dbSNP,vcf,exact,0,KRG,K1,dbSNPm \
# --custom ${Resource}clinvar.vcf.gz,Clinvar.vcf.gz,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
# --plugin dbNSFP,/home/goldpm1/tools/ANNOVAR/humandb/dbnsfp4.2a/dbNSFP4.0a.gz,ALL \
# --plugin SpliceAI,snv=${Resource}spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=${Resource}spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
# --cache \
# --dir_cache /data/public/VEP/104 \
# --cache_version 104 \
# --vcf \
# --offline



/opt/Yonsei/ensembl-vep/101.0/vep \
-v -assembly "GRCh38" --everything --terms "SO" --fork "24" \
-i ${INPUT_VCF} \
--format "vcf" \
-o ${OUTPUT_VCF} \
--force --no_stats \
--vcf \
--dir_plugins /data/public/VEP/101/Plugins \
--custom ${Resource}dbSNP/GCF_000001405.39.re.vep_custom2.vcf.gz,dbSNP,vcf,exact,0,KRG,K1,dbSNPm \
--custom ${Resource}clinvar.vcf.gz,Clinvar.vcf.gz,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
--fasta $REF \
--cache_version 101 \
--dir_cache /data/public/VEP/101 \
--plugin SpliceAI,snv=/data/project/DC_WGS/Resource/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/data/project/DC_WGS/Resource/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
--plugin dbNSFP,/home/goldpm1/tools/ANNOVAR/humandb/dbnsfp4.2a/dbNSFP4.0a.gz,ALL \
--plugin AlphaMissense,file=/data/resource/annotation/human/Alphamissense/AlphaMissense_hg38.tsv.gz \
--offline \
--cache



# ############ NEAREST GENE annotation #####################################

# # 일단 header만 따로 떼낸다
# echo -e "python3 /data/project/craniosynostosis/script/noncoding/noncoding_pipe_6_addheader.py --INPUTVCF ${OUTPUT_VCF} --OUTPUT_VCF ${OUTPUT_VCF}header"
# python3 "/data/project/craniosynostosis/script/noncoding/noncoding_pipe_6_addheader.py" --INPUTVCF ${OUTPUT_VCF} --OUTPUT_VCF ${OUTPUT_VCF}"header"

# # bedtools closest 시행. 뒤에 5개가 달라붙는다
# GENE_BED="/home/goldpm1/resources/whole.genelist.hg38.bedtools_merge.bed"
# echo -e "bedtools closest -a ${OUTPUT_VCF}  -b ${GENE_BED} > ${OUTPUT_VCF}temp"
# bedtools closest -a ${OUTPUT_VCF}  -b ${GENE_BED} > ${OUTPUT_VCF}"temp"

# # 다시 parsing 해서 집어넣는다
# echo -e "python3 /data/project/craniosynostosis/script/noncoding/noncoding_pipe_6_vcfparsing.py --INPUTVCF ${OUTPUT_VCF}temp --OUTPUT_VCF ${OUTPUT_VCF}body"
# python3 "/data/project/craniosynostosis/script/noncoding/noncoding_pipe_6_vcfparsing.py" --INPUTVCF ${OUTPUT_VCF}"temp" --OUTPUT_VCF ${OUTPUT_VCF}"body"
# cat ${OUTPUT_VCF}"header" ${OUTPUT_VCF}"body" > ${OUTPUT_VCF}
# rm -rf  ${OUTPUT_VCF}"header" ${OUTPUT_VCF}"body" ${OUTPUT_VCF}"temp"

bgzip -c -f ${OUTPUT_VCF} > ${OUTPUT_VCF}".gz"
tabix -p vcf ${OUTPUT_VCF}".gz"
