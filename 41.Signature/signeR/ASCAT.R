library ("ASCAT")
library ("GenomicRanges")

ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", 
                          Tumor_BAF_file = "Tumor_BAF.txt", 
                          Germline_LogR_file = "Germline_LogR.txt",
                          Germline_BAF_file = "Germline_BAF.txt", 
                          gender = rep('XX',100), 
                          genomeVersion = "hg38")

ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")

ascat.bc = ascat.correctLogR(ascat.bc, 
                             GCcontentfile = "GC_example.txt", 
                             replictimingfile = "RT_example.txt")

ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")


ascat.bc = ascat.asmultipcf(ascat.bc)

ascat.plotSegmentedData(ascat.bc)

ascat.output = ascat.runAscat(ascat.bc, 
                              write_segments = T) # gamma=1 for HTS data

QC = ascat.metrics(ascat.bc,ascat.output)

save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')