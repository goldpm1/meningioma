library("MutationTimeR")
library("GenomicRanges")
library ("stringr")
library("argparse")

parser <- ArgumentParser()

parser$add_argument("--Sample_ID")
parser$add_argument("--TISSUE")
parser$add_argument("--SEQUENZA_SEGMENT_PATH")
parser$add_argument("--SEQUENZA_PLOIDY_PATH")
parser$add_argument("--MUTATIONTIMER_INPUT_VCF_PATH")
parser$add_argument("--MUTATIONTIMER_RESULT_DIR")

args <- parser$parse_args()

Sample_ID  = args$Sample_ID 
TISSUE  = args$TISSUE 
SEQUENZA_SEGMENT_PATH  = args$SEQUENZA_SEGMENT_PATH 
SEQUENZA_PLOIDY_PATH = args$SEQUENZA_PLOIDY_PATH
MUTATIONTIMER_INPUT_VCF_PATH = args$MUTATIONTIMER_INPUT_VCF_PATH
MUTATIONTIMER_RESULT_DIR  = args$MUTATIONTIMER_RESULT_DIR 


paste0 ("SEQUENZA_PLOIDY_PATH : ", SEQUENZA_PLOIDY_PATH)
paste0 ("MUTATIONTIMER_INPUT_VCF_PATH : ", MUTATIONTIMER_INPUT_VCF_PATH)
paste0 ("MUTATIONTIMER_RESULT_DIR  : ", MUTATIONTIMER_RESULT_DIR )


ploidy_df = read.table ( SEQUENZA_PLOIDY_PATH, header = T, sep = "\t")
TUMOR_PURITY = ploidy_df[1,"cellularity"]


# Mutation 정보
meningioma_vcf = VariantAnnotation::readVcf( MUTATIONTIMER_INPUT_VCF_PATH )


# CNV 정보
df = read.table (SEQUENZA_SEGMENT_PATH, header = T, sep = "\t")
df$chromosome <- sub("^chr", "", df$chromosome)  # 맨 앞에 달려있는 chr을 빼주기
gr = GRanges(seqnames = df$chromosome, 
             ranges = IRanges (df$start.pos, df$end.pos), 
             strand = rep("*", nrow(df)),
             major_cn = df$A, minor_cn = df$B, 
             clonal_frequency = TUMOR_PURITY )



clusters = data.frame(cluster = 2, proportion=TUMOR_PURITY, n_ssms=100)
mt = MutationTimeR::mutationTime(meningioma_vcf, gr, n.boot=10)


meningioma_vcf = addMutTime(meningioma_vcf, mt$V)
writeVcf(meningioma_vcf, paste0( MUTATIONTIMER_RESULT_DIR, "/", Sample_ID, "_", TISSUE, ".MutationTimeR.vcf") )
mcols(gr) = cbind(mcols(gr), mt$T)

options(bitmapType="cairo")
png( paste0(MUTATIONTIMER_RESULT_DIR, "/", Sample_ID, "_", TISSUE, ".results.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 300 )
plotSample (vcf = meningioma_vcf, cn = gr, title = paste0 (Sample_ID, "_", TISSUE ,  "  (tumor purity = ", TUMOR_PURITY, ")" )  )
dev.off()