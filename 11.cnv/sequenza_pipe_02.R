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


print('Sequenza.extract....')
sequenza_extract = sequenza.extract(SEQUENZA_SMALL_SEQZ, verbose =TRUE, normalization.method="median", 
                                                                assembly = hg, chromosome.list = chrom)

# sequenza_extract = sequenza.extract("/home/goldpm1/Meningioma/11.cnv/1.seqz/220930_Tumor.small.seqz.gz", normalization.method="median", assembly='hg38')
                                                            


print( 'Sequenza.fit....')
CP = sequenza.fit (sequenza_extract, 
                                ploidy = seq(1.5, 2.5, 0.25), 
                                cellularity = seq(0.3,1,0.2),
                                )

print( 'Sequenza.results...')
sequenza.results(sequenza.extract = sequenza_extract, cp.table = CP, sample.id = ID, out.dir= SEQUENZA_OUTPUT_DIR, CNt.max = 5 )   




# print('Sequenza.extract....')
# sequenza_extract = sequenza.extract(SEQUENZA_SMALL_SEQZ, verbose =TRUE, normalization.method="median", 
#                                                                 assembly='hg38', chromosome.list = chrom, 
#                                                                 min.reads = 10, min.reads.normal = 10, min.reads.baf = 10, 
#                                                                 mufreq.treshold = 0.15)

# print( 'Sequenza.fit....')
# CP = sequenza.fit (sequenza_extract, 
#                                 ploidy = seq(1.5, 2.5, 0.25), 
#                                 cellularity = seq(0.3,1,0.2),
#                                 segment.filter = 3e3, 