#!/bin/bash
#$ -S /bin/bash
#$ -cwd

### 실행code:  bash /data/project/Meningioma/script/03.Variant_calling\&Annotation/Multiple_pipe_04.Sequoia.sh

module load mpboot

# "190426", "230323", "241127"
Sample_ID="190426"
OUTPUT_DIR="/data/project/Meningioma/04.mutect/09.multiple/"${Sample_ID}"/Sequoia"
if [ ! -d ${OUTPUT_DIR} ] ; then
    mkdir ${OUTPUT_DIR}
fi


#1. Sequoia_mpboot 돌리기 (여기서는 germline을 다 제거했으니까 cutoff를 극단적으로 높여도 됨)
Rscript /home/goldpm1/tools/Sequoia/build_phylogeny.R \
-r "/data/project/Meningioma/04.mutect/09.multiple/"${Sample_ID}"/depth.tsv" \
-v "/data/project/Meningioma/04.mutect/09.multiple/"${Sample_ID}"/alt.tsv" \
-o ${OUTPUT_DIR}"/" \
--germline_cutoff -1 \
-b TRUE \
--only_snvs TRUE \
-i ${Sample_ID}

#--snv_rho 0.0000001 \

#2. ggtree로 그리기
for OPTION in "both" "snv" "indel"; do
    Rscript "/data/project/Meningioma/script/61.Lowinput/08.Sequoia_pipe_03.ggtree.R" \
        ${OUTPUT_DIR}"/"${Sample_ID}"_"${OPTION}"_tree_with_branch_length.tree" \
        ${OUTPUT_DIR}"/"${Sample_ID}"_ggtree_"${OPTION}".pdf"
done