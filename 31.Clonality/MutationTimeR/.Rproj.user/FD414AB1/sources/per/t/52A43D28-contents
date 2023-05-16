library("MutationTimeR") 
library ("BSgenome.Hsapiens.UCSC.hg38")
DIR = "/home/goldpm1/Meningioma/31.Clonality/03.mutationtimeR/SeYoung"

SY_vcf <- readVcf(paste0(DIR, "/P26_tumor.RGadded.marked.fixed.mutect2.filter.PASS.ann.timer.vcf")) 

Sequenza.segments <- read.csv(paste0( DIR, "/P26_tumor.RGadded.marked.fixed.small.seqz_segments.txt"), sep = '\t', header = TRUE)
Sequenza.alernative_solutions <- read.csv(paste0(DIR, "/P26_tumor.RGadded.marked.fixed.small.seqz_alternative_solutions.txt"), sep = '\t', header = TRUE)
purity = 1 - Sequenza.alernative_solutions[1,1]

SY_bb <- GRanges(seqnames = Sequenza.segments$chromosome,
              ranges = IRanges(start = Sequenza.segments$start.pos,
                               end = Sequenza.segments$end.pos),
              major_cn = Sequenza.segments$A,
              minor_cn = Sequenza.segments$B,
              clonal_frequency = purity) # Copy number segments, needs columns  major_cn, minor_cn and clonal_frequency of each segment

clusters = data.frame(cluster = 1,
                      proportion=purity,
                      n_ssms=100)

SY_mt <- mutationTime(SY_vcf, SY_bb, clusters, n.boot=10)

SY_vcf <- addMutTime(vcf, SY_mt$V)

mcols(SY_bb) <- cbind(mcols(SY_bb),mt$T)
#export.bed(bb, file = paste0("TitanCNA/50kb/",i,"_tumor.segs.Mut.Time.bed"))

pdf(paste0("result_Sequenza/",i,"_ventricle_B.MutationTimeR.pdf"), width = 10, height = 8)
plotSample(SY_vcf, SY_bb)
dev.off()

SY_vcf@metadata
