library("MutationTimeR")
library("GenomicRanges")
library ("stringr")
library("argparse")

parser <- ArgumentParser()

parser$add_argument("--Sample_ID")
parser$add_argument("--TISSUE")
parser$add_argument("--FACET_CNV_MATRIX_PATH")
parser$add_argument("--MUTATIONTIMER_INPUT_VCF_PATH")
parser$add_argument("--MUTATIONTIMER_RESULT_DIR")

args <- parser$parse_args()

Sample_ID  = args$Sample_ID 
TISSUE  = args$TISSUE 
FACET_CNV_MATRIX_PATH  = args$FACET_CNV_MATRIX_PATH 
MUTATIONTIMER_INPUT_VCF_PATH = args$MUTATIONTIMER_INPUT_VCF_PATH
MUTATIONTIMER_RESULT_DIR  = args$MUTATIONTIMER_RESULT_DIR 


paste0 ("FACET_CNV_MATRIX_PATH : ", FACET_CNV_MATRIX_PATH)
paste0 ("MUTATIONTIMER_INPUT_VCF_PATH : ", MUTATIONTIMER_INPUT_VCF_PATH)
paste0 ("MUTATIONTIMER_RESULT_DIR  : ", MUTATIONTIMER_RESULT_DIR )


facetcnv_df = read.table ( FACET_CNV_MATRIX_PATH, header = T, sep = "\t")
TUMOR_PURITY = facetcnv_df[1,"TUMOR_PURITY"]
# Create a GRange object  (chr을 빼주자)
gr = GRanges(seqnames = str_replace (facetcnv_df$CHR, pattern = "chr", replacement = "" ), 
             ranges = IRanges (facetcnv_df$START, facetcnv_df$END), 
             strand = rep("*", nrow(facetcnv_df)),
             major_cn = facetcnv_df$MAJOR_CN, minor_cn = facetcnv_df$MINOR_CN, 
             clonal_frequency = TUMOR_PURITY )

meningioma_vcf = VariantAnnotation::readVcf( MUTATIONTIMER_INPUT_VCF_PATH )
clusters = data.frame(cluster = 1, proportion=TUMOR_PURITY,n_ssms=100)
mt = MutationTimeR::mutationTime(meningioma_vcf, gr, n.boot=10)


meningioma_vcf <- addMutTime(meningioma_vcf, mt$V)
writeVcf(meningioma_vcf, paste0( MUTATIONTIMER_RESULT_DIR, "/", Sample_ID, "_", TISSUE, ".MutationTimeR.vcf") )
mcols(gr) <- cbind(mcols(gr),mt$T)

options(bitmapType="cairo")
png( paste0(MUTATIONTIMER_RESULT_DIR, "/", Sample_ID, "_", TISSUE, ".results.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 300 )
plotSample (vcf = meningioma_vcf, cn = gr, title = paste0 (Sample_ID, "_", TISSUE, "  (tumor purity = ", TUMOR_PURITY, ")") )
dev.off()