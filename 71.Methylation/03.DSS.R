library(DSS)

DIR="/data/project/Meningioma/71.Methylation/03.call"

dat1.1 = read.table(file.path(DIR, "230526_Tumor", "230526_Tumor.DSS.txt"), header=TRUE)
dat1.2 = read.table(file.path(DIR, "221026_Tumor", "221026_Tumor.DSS.txt"), header=TRUE)
#dat1.3 = read.table(file.path(DIR, "221102_Tumor", "221102_Tumor.DSS.txt"), header=TRUE)
dat1.4 = read.table(file.path(DIR, "230323_Tumor", "230323_Tumor.DSS.txt"), header=TRUE)

BSobj <- makeBSseqData( list( dat1.1, dat1.2, dat1.4), c("230526", "221026" , "230323"))

BSobj <- BSmooth(BSobj, ns=70, h=1000, maxGap=100000)



dmlTest = DMLtest(BSobj, group1=c("230526"), group2=c("221026", "230323"))  # There is no biological replicates in at least one condition. Please set smoothing=TRUE or equal.disp=TRUE and retry.
dmlTest.sm = DMLtest(BSobj, group1=c("230526"), group2=c("221026", "230323"), smoothing=TRUE)

# Find DML (Loci) & DMR (Region)
dmls <- callDML(dmlTest.sm, p.threshold=0.01)
dmrs <- callDMR(dmlTest.sm, p.threshold=0.01, minCG=3, dis.merge=100)

showOneDMR(dmrs[1,], BSobj)

write.table(dmrs, file="DMRs_DSS.txt", sep="\t", row.names=FALSE, col.names=TRUE)

dmrs_chr22 = dmrs[dmrs$chr=="chr22",]
View ( dmrs_chr22 [ order (dmrs_chr22$start), ] )
