#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load mpboot

OUTPUT_DIR="/data/project/Meningioma/04.mutect/07.2D_merged/04.Sequoia"
if [ ! -d ${OUTPUT_DIR} ] ; then
    mkdir ${OUTPUT_DIR}
fi

Rscript /home/goldpm1/tools/Sequoia/build_phylogeny.R \
-r "/data/project/Meningioma/04.mutect/07.2D_merged/02.BCFTOOLS_MERGE_VCF/190426.merged.heatmap.depth.tsv" \
-v "/data/project/Meningioma/04.mutect/07.2D_merged/02.BCFTOOLS_MERGE_VCF/190426.merged.heatmap.alt.tsv" \
-o "/data/project/Meningioma/04.mutect/07.2D_merged/04.Sequoia/" \
--germline_cutoff -30 \
--only_snvs TRUE \
-i "190426"


# Rscript "08.Sequoia_ggtree.R" \
# "/data/project/Meningioma/61.Lowinput/01.XT_HS/06.HC/07.2D_merged/03.Sequoia/190426_both_tree_with_branch_length.tree" \
# "/data/project/Meningioma/61.Lowinput/01.XT_HS/06.HC/07.2D_merged/03.Sequoia/190426_both_tree_ggtree.jpg"