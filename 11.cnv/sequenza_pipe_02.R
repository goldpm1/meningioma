library(argparse)
library(sequenza)
library(copynumber)

# Input by argparser
parser <- ArgumentParser()
parser$add_argument("--ID", default = "")
parser$add_argument("--hg", default = "")
parser$add_argument("--SEQUENZA_SMALL_SEQZ", default = "")
parser$add_argument("--SEQUENZA_OUTPUT_DIR", default = "")
args <- parser$parse_args()

ID = args$ID
hg = args$hg
SEQUENZA_OUTPUT_DIR = args$SEQUENZA_OUTPUT_DIR
SEQUENZA_SMALL_SEQZ = args$SEQUENZA_SMALL_SEQZ


print (ID)
print (hg)
print (SEQUENZA_SMALL_SEQZ)
print (SEQUENZA_OUTPUT_DIR)


args=commandArgs(trailingOnly=TRUE)


chrom = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
chrom = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")


print('Sequenza.extract....')  # 원래는 weighted.mean = TRUE, window = 1e7
sequenza_extract = sequenza.extract( SEQUENZA_SMALL_SEQZ, verbose =TRUE, normalization.method="median", weighted.mean = FALSE, window = 6e7, 
                                                                assembly = hg, chromosome.list = chrom)
                                                    


print( 'Sequenza.fit....')    # segment.filter = 3e6, 원래는 N.ratio.filter = 10,
CP = sequenza.fit (sequenza_extract, 
                                segment.filter = 6e7, N.ratio.filter = 50,
                                ploidy = seq(2.0, 2.5, 0.25), 
                                cellularity = seq(0.3,1,0.2)
                                )

print( 'Sequenza.results...')     # CNt.max   total allele count max가 3인걸 뽑느다
sequenza.results(sequenza.extract = sequenza_extract, cp.table = CP, sample.id = ID, out.dir= SEQUENZA_OUTPUT_DIR, CNt.max = 3 ) 