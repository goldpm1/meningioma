library("signeR")
library("VariantAnnotation")
library("BSgenome.Hsapiens.UCSC.hg38")
library("rtracklayer")

CURRENT_PATH=getwd()

vcf_files = Sys.glob("/home/goldpm1/Meningioma/04.mutect/05.unique/*.vcf")

#01. Make sample*96 matrix from multi-sample vcfs
mut = matrix(ncol=96,nrow=0)
for(i in vcf_files) {
  vo = readVcf(i, "hg38")
  # colnames(vo) = i        # 이렇게 안하면 자동으로 prefix가 된다
  m0 = signeR::genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg38, vo)
  mut = rbind(mut, m0)
}


#02. make opportunity matrix
target_regions <- rtracklayer::import(con="/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/a.bed", format="bed")
opp = signeR::genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg38, target_regions, nsamples=nrow(mut))


#03. Load reference COSMIC database
Pmatrix_BRC = as.matrix(read.table(system.file("extdata","Cosmic_signatures_BRC.txt", package="signeR"), sep="\t", check.names=FALSE))      # 96 context * 6 signatures
Pmatrix_COSMICv3 = as.matrix(read.table( "/home/goldpm1/resources/COSMIC_SBS/COSMIC_v3.3.1_SBS_GRCh38.txt", sep="\t", header = TRUE, check.names=FALSE, row.names = 1))   # 96 context * 80 signatures
#rownames()


#04.signature decomposition
#signatures <- signeR::signeR(M=mut, Opport=opp,  P=Pmatrix_COSMICv3 , fixedP=FALSE, main_eval=100, EM_eval=50, EMit_lim=20, nlim=c(2,4))
signatures <- signeR::signeR(M=mut, Opport=opp, EMit_lim=20, nlim=c(2,4))

#05. 최적의 K값을 BIC plot을 통해 찾기
options(bitmapType="cairo")
png( paste0(CURRENT_PATH, "/BICboxplot.png") )
signeR::BICboxplot (signatures)
dev.off()

#06. Results and Plots
options(bitmapType="cairo")
png( paste0(CURRENT_PATH, "/path.png") )
signeR::Paths(signatures$SignExposures)
dev.off()

options(bitmapType="cairo")
png( paste0(CURRENT_PATH, "/SignPlot.png") )
signeR::SignPlot(signatures$SignExposures)
dev.off()

options(bitmapType="cairo")
png( paste0(CURRENT_PATH, "/SignHeat.png") )
signeR::SignHeat(signatures$SignExposures)
dev.off()

options(bitmapType="cairo")
png( paste0(CURRENT_PATH, "/ExposureBoxplot.png") )
signeR::ExposureBoxplot(signatures$SignExposures)
dev.off()

options(bitmapType="cairo")
png( paste0(CURRENT_PATH, "/ExposureBarplot.png") )
signeR::ExposureBarplot(signatures$SignExposures, relative=TRUE)
dev.off()

options(bitmapType="cairo")
png( paste0(CURRENT_PATH, "/ExposureHeat.png") )
signeR::ExposureHeat(signatures$SignExposures)
dev.off()
