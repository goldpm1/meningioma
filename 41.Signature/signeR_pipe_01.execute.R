library("signeR")
packageVersion("signeR")
library("VariantAnnotation")
library("BSgenome.Hsapiens.UCSC.hg38")
library("rtracklayer")
library("argparse")

parser <- ArgumentParser()

parser$add_argument("--SIGNER_INPUT_VCF_DIR")
parser$add_argument("--SIGNER_RESULT_DIR")
parser$add_argument("--SIGNER_RESULT_Phat_PATH")

args <- parser$parse_args()

SIGNER_INPUT_VCF_DIR = args$SIGNER_INPUT_VCF_DIR
SIGNER_RESULT_DIR = args$SIGNER_RESULT_DIR 
SIGNER_RESULT_Phat_PATH = args$SIGNER_RESULT_Phat_PATH

CURRENT_PATH=getwd()
PARENT_PATH=paste0( paste0( unlist(stringr::str_split(getwd(), "/")) [1: length ( unlist(stringr::str_split(getwd(), "/"))) - 1], collapse = "/" ) )
RESULT_PATH=paste0( SIGNER_RESULT_DIR )


paste0 ("SIGNER_INPUT_VCF_DIR : ", SIGNER_INPUT_VCF_DIR)
paste0 ("SIGNER_RESULT_DIR : ", SIGNER_RESULT_DIR)





vcf_files = Sys.glob( paste0( SIGNER_INPUT_VCF_DIR, "/*.vcf")  )

#01. Make sample*96 matrix from multi-sample vcfs
mut = matrix(ncol=96,nrow=0)
for(i in vcf_files) {
  vo = readVcf(i, "hg38")
  # colnames(vo) = i        # 이렇게 안하면 자동으로 prefix가 된다
  m0 = signeR::genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg38, vo)
  mut = rbind(mut, m0)
}


#02. make opportunity matrix
target_regions <- rtracklayer::import(con= paste0("/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/a.bed"), format="bed")
opp = signeR::genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg38, target_regions, nsamples=nrow(mut))


#03. Load reference COSMIC database
#Pmatrix_BRC = as.matrix(read.table(system.file("extdata","Cosmic_signatures_BRC.txt", package="signeR"), sep="\t", check.names=FALSE))      # 96 context * 6 signatures
Pmatrix_COSMICv3 = as.matrix(read.table( paste0("/home/goldpm1/resources/COSMIC_SBS/COSMIC_v3.3.1_SBS_GRCh38.txt"), sep="\t", header = TRUE, check.names=FALSE, row.names = 1))   # 96 context * 80 signatures


#04.signature decomposition
#signatures <- signeR::signeR(M=mut, Opport=opp, P=Pmatrix_COSMICv3 , fixedP=FALSE, main_eval=100, EM_eval=50, EMit_lim=20, nlim=c(2,4))
signatures <- signeR::signeR(M=mut, Opport=opp, EMit_lim=20, nlim=c(2,4))





#05. 최적의 K값을 BIC plot을 통해 찾기
options(bitmapType="cairo")
png( paste0(RESULT_PATH, "/BICboxplot.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 300 )
signeR::BICboxplot (signatures)
dev.off()

#06. Results and Plots
write.table ( as.data.frame(signatures$Phat) ,  file = SIGNER_RESULT_Phat_PATH, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

options(bitmapType="cairo")
png( paste0(RESULT_PATH, "/path.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 100 )
signeR::Paths(signatures$SignExposures)
dev.off()

options(bitmapType="cairo")
png( paste0(RESULT_PATH, "/SignPlot.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 100 )
signeR::SignPlot(signatures$SignExposures)
dev.off()

options(bitmapType="cairo")
png( paste0(RESULT_PATH, "/SignHeat.png"), width=1080, height=1900, units = "px", pointsize = 10, res = 100 )
signeR::SignHeat(signatures$SignExposures)
dev.off()

options(bitmapType="cairo")
png( paste0(RESULT_PATH, "/ExposureBoxplot.png"), width=1920, height=1080, units = "px", pointsize = 8, res = 100 )
signeR::ExposureBoxplot(signatures$SignExposures)
dev.off()

options(bitmapType="cairo")
png( paste0(RESULT_PATH, "/ExposureBarplot_absolute.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 100 )
signeR::ExposureBarplot(signatures$SignExposures, relative=FALSE)
dev.off()

options(bitmapType="cairo")
png( paste0(RESULT_PATH, "/ExposureBarplot_relative.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 100 )
signeR::ExposureBarplot(signatures$SignExposures, relative=TRUE)
dev.off()

options(bitmapType="cairo")
png( paste0(RESULT_PATH, "/ExposureHeat.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 100 )
signeR::ExposureHeat(signatures$SignExposures)
dev.off()

png( paste0(RESULT_PATH, "/Hieraquical_Clustering.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 100 )
HCE <- HClustExp(signatures$SignExposures,method.dist="euclidean", method.hclust="average")
dev.off()

png( paste0(RESULT_PATH, "/Fuzzy_Clustering.png"), width=1920, height=1080, units = "px", pointsize = 10, res = 100 )
FCE <- FuzzyClustExp(signatures$SignExposures, Clim=c(3,7))
dev.off()