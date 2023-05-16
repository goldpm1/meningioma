library("MutationTimeR")
library("GenomicRanges")
library ("stringr")

packageVersion("signeR")

CURRENT_PATH=getwd()

# INPUT_DF
SEQUENZA_SEGMENT_PATH="/home/goldpm1/Meningioma/11.cnv/2.sequenza/221026_Tumor_segments.txt"
SEQUENZA_PLOIDY_PATH="/home/goldpm1/Meningioma/11.cnv/2.sequenza/221026_Tumor_confints_CP.txt"
MUTECT_RESCUE_PATH="/home/goldpm1/Meningioma/31.Clonality/03.mutationtimeR/01.vcf/221026_Tumor.MT2.FMC.HF.RMBLACK.vep.rescue.removechr.vcf"
HC_PATH="/home/goldpm1/Meningioma/06.hc/03.HF/221026/Tumor/221026_Tumor.DP100.vcf"


ploidy_df = read.table ( SEQUENZA_PLOIDY_PATH, header = T, sep = "\t")
TUMOR_PURITY = ploidy_df[1,"cellularity"]

df = read.table (SEQUENZA_SEGMENT_PATH, header = T, sep = "\t")


# Create a GRange object  (chr을 빼주자)
gr = GRanges(seqnames = str_replace (df$chromosome, pattern = "chr", replacement = "" ), 
             ranges = IRanges (df$start.pos, df$end.pos), 
             strand = rep("*", nrow(df)),
             major_cn = df$A, minor_cn = df$B, 
             clonal_frequency = TUMOR_PURITY )

meningioma_vcf = VariantAnnotation::readVcf( MUTECT_RESCUE_PATH )
clusters = data.frame(cluster = 1, proportion=TUMOR_PURITY,n_ssms=100)
mt = MutationTimeR::mutationTime(meningioma_vcf, gr, n.boot=10)


meningioma_vcf <- addMutTime(meningioma_vcf, mt$V)
writeVcf(meningioma_vcf, paste0( CURRENT_PATH, "/YS_mningioma.vcf") )
mcols(gr) <- cbind(mcols(gr),mt$T)


options(bitmapType="cairo")
png( paste0("/home/goldpm1/Meningioma/script/31.Clonality/MutationTimeR/YS_results.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 300 )
plotSample (vcf = meningioma_vcf, cn = gr, title = "221026_Tumor")
dev.off()




###############################################################################################



FACET_CNV_MATRIX_PATH="/home/goldpm1/Meningioma/31.Clonality/01.make_matrix/221026/221026_Tumor.facetcnv_to_bed_df.tsv"
MUTATIONTIMER_INPUT_VCF_PATH = "/home/goldpm1/Meningioma/31.Clonality/03.mutationtimeR/01.vcf/221026_Tumor.MT2.FMC.HF.RMBLACK.vep.rescue.removechr.vcf"
facetcnv_df = read.table ( FACET_CNV_MATRIX_PATH, header = T, sep = "\t")
TUMOR_PURITY = facetcnv_df[1,"TUMOR_PURITY"]
gr = GRanges(seqnames = str_replace (facetcnv_df$CHR, pattern = "chr", replacement = "" ), 
             ranges = IRanges (facetcnv_df$START, facetcnv_df$END), 
             strand = rep("*", nrow(facetcnv_df)),
             major_cn = facetcnv_df$MAJOR_CN, minor_cn = facetcnv_df$MINOR_CN, 
             clonal_frequency = TUMOR_PURITY )
facetcnv_df = read.table ( FACET_CNV_MATRIX_PATH, header = T, sep = "\t")

meningioma_vcf = VariantAnnotation::readVcf( MUTATIONTIMER_INPUT_VCF_PATH )
clusters = data.frame(cluster = 3, proportion=TUMOR_PURITY,n_ssms=100)
mt = MutationTimeR::mutationTime(meningioma_vcf, gr, n.boot=10)


meningioma_vcf <- addMutTime(meningioma_vcf, mt$V)
writeVcf(meningioma_vcf, paste0( MUTATIONTIMER_RESULT_DIR, "/", Sample_ID, "_", TISSUE, ".MutationTimeR.vcf") )
mcols(gr) <- cbind(mcols(gr),mt$T)

options(bitmapType="cairo")
png( paste0(MUTATIONTIMER_RESULT_DIR, "/", Sample_ID, "_", TISSUE, ".results.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 300 )
plotSample (vcf = meningioma_vcf, cn = gr, title = paste0 (Sample_ID, "_", TISSUE, "  (tumor purity = ", TUMOR_PURITY, ")") )
dev.off()


