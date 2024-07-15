#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load mpboot

# 190426
Rscript /home/goldpm1/tools/Sequoia/build_phylogeny.R \
-r "/data/project/Meningioma/61.Lowinput/01.XT_HS/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/190426.heatmap.depth.tsv" \
-v "/data/project/Meningioma/61.Lowinput/01.XT_HS/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/190426.heatmap.alt.tsv" \
-o "/data/project/Meningioma/61.Lowinput/01.XT_HS/06.HC/07.2D_merged/03.Sequoia/" \
--germline_cutoff -30 \
--only_snvs TRUE \
-i "190426"


Rscript "08.Sequoia_ggtree.R" \
"/data/project/Meningioma/61.Lowinput/01.XT_HS/06.HC/07.2D_merged/03.Sequoia/190426_both_tree_with_branch_length.tree" \
"/data/project/Meningioma/61.Lowinput/01.XT_HS/06.HC/07.2D_merged/03.Sequoia/190426_both_tree_ggtree.jpg"



# # 230405
Rscript /home/goldpm1/tools/Sequoia/build_phylogeny.R \
-r "/data/project/Meningioma/61.Lowinput/02.PTA/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/230405.heatmap.depth.tsv" \
-v "/data/project/Meningioma/61.Lowinput/02.PTA/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/230405.heatmap.alt.tsv" \
-o "/data/project/Meningioma/61.Lowinput/02.PTA/06.HC/07.2D_merged/03.Sequoia/" \
--germline_cutoff -30 \
--only_snvs TRUE \
-i "230405"


Rscript "08.Sequoia_ggtree.R" \
"/data/project/Meningioma/61.Lowinput/02.PTA/06.HC/07.2D_merged/03.Sequoia/230405_both_tree_with_branch_length.tree" \
"/data/project/Meningioma/61.Lowinput/02.PTA/06.HC/07.2D_merged/03.Sequoia/230405_both_tree_ggtree.jpg"
