import pysam, glob
import pandas as pd
import numpy as np

class Variant:
    def __init__ (self, date, bamtype, BAM_Dir, chr_, pos_):
        self.date = date
        self.bamtype = bamtype
        self.BAM_Dir =  BAM_Dir
        self.chr_ = str ( chr_ )
        self.pos_ = int (pos_)
    def pysam_read (self):
        samfile = pysam.AlignmentFile(self.BAM_Dir, 'rb')
        #for ID in _variant_dic:
        multiallele_check = set()
        ins_count, ins_list, del_count, del_list = 0, [], 0, []
        base_counts = {"A": 0, "T": 0, "C": 0, "G": 0, "N" : 0, "Ins" : 0 , "Del": 0}

        iter = samfile.fetch(contig = self.chr_, start = self.pos_-1, end = self.pos_)
        for i, read in enumerate(iter):
            if read.is_proper_pair:
                if read.flag > 512:
                    continue
                read_info = str(read).rstrip().split('\t')
                cigar_string = read.cigarstring
                start_site = read.reference_start + 1
                end_site = read.reference_end
                paired_read_start_site = read.next_reference_start + 1
                read_base = read.query_sequence
                read_length = read.query_alignment_length
                soft_al = read.query_alignment_start ### without clipped base
                soft_al2 = read.query_alignment_end ### without clipped base

                if ('H' in cigar_string) | ('S' in cigar_string):
                    continue
                if (self.pos_ >= read.reference_end - 4) | (self.pos_ <= read.reference_start + 4):
                    continue
            
                if 'D' in cigar_string:
                    if start_site > self.pos_ :
                        print (cigar_string, self.pos_, start_site)
                    pp = start_site
                    for j in cigar_string.split("D")[0].split("M")[0:-1]:
                        if ("I" in j) | ("S" in j):
                            continue
                        else:
                            pp += int (j)
                    if pp - 1== int(self.pos_):
                        base_counts["Del"] += 1
                        #del_list.append (read.query_name)
                        del_list.append (read.cigarstring)
                elif 'I' in cigar_string:
                    if start_site > self.pos_ :
                        print (cigar_string, self.pos_, start_site)
                    pp = start_site
                    for j in cigar_string.split("I")[0].split("M")[0:-1]:
                        if ("D" in j) | ("S" in j):
                            continue
                        else:
                            pp += int (j)
                    if pp - 1== int(self.pos_):
                        base_counts["Ins"] += 1
                        #ins_list.append (read.query_name)
                        ins_list.append (read.cigarstring)
                else:
                    base_quality = read.query_qualities[ self.pos_ - start_site ]
                    if base_quality >= 15:   #12로 하면 230920을 살리고, 15로 하면 230419를 좋게 만들고
                        base_counts [ read.query_alignment_sequence [ self.pos_ - start_site ] ] += 1
                        # if base_counts [ read.query_alignment_sequence [ self.pos_ - start_site ] ] == 1:
                        #     print (  start_site,  cigar_string, end = "\t"  )
                        #     for i in range ( self.pos_- 2, self.pos_ + 3 ):
                        #         print (  read.query_alignment_sequence [ i - start_site ] , end = "" ) 
                        #     print ( "" )
                            
                    continue

        print ("\t# of total read depth : {}". format ( np.sum ( list ( base_counts.values() )  ) ) )
        for key in base_counts.keys():
            print ("\t# of {} : {}". format (key, base_counts[key]))
        print ("\t\tins_list = {}".format (ins_list))
        print ("\t\tdel_list = {}".format (del_list))




interval_df = pd.read_csv ("/data/project/Meningioma/00.Targetted/07.pysam/interval.bed", header = None, sep = "\t")
samplename_list =  list ( pd.read_csv ("/data/project/Meningioma/00.Targetted/07.pysam/sample_name.txt", header = None, sep = "\t").iloc[:, 0] )
BAM_PATH_LIST = sorted ( glob.glob ( "/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/01.Pre_bam/*.bam" ) )

interval_df ["ID"]  = interval_df.iloc[:, 0].astype(str) + ":" + interval_df.iloc[:, 1].astype(str) + "(" + interval_df.iloc[:, 4].astype(str) + ")"

for DATE in samplename_list:
    SAMPLENAME = DATE
    BAM_PATH = "/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/01.Pre_bam/" + DATE + ".sorted.bam"
    #BAM_PATH = "/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/06.Final_bam/" + DATE + "_Dura.bam"

    count_list = []

    for k in range (0,  interval_df.shape[0]):
        CHR, POS, REF, ALT, GENE, TYPE = interval_df.iloc[k][0], interval_df.iloc[k][1], interval_df.iloc[k][2], interval_df.iloc[k][3], interval_df.iloc[k][4], interval_df.iloc[k][5]

        ALT_TYPE = TYPE if TYPE in ["Ins", "Del"] else ALT

        print ( "\n{} ( {}:{} {}>{} ) ".format ( DATE, CHR, POS, REF, ALT))
        Count = Variant ( DATE, "Amplicon", BAM_PATH, CHR, POS  )
        Count.pysam_read()




#BAM_Dir, chr_, pos_ = "/data/project/Meningioma/02.Align/hg38/Dura/05.Final_bam/230419_Dura_no_softclips.bam", "chr22", 29674838                  # WES




# Dura_230323_WES1 = Variant ( 230323, "WES1", "/data/project/Meningioma/02.Align/hg38/Dura/05.Final_bam/230323_Dura.bam", "chr22", 29661335  )
# Dura_230323_WES2 = Variant ( 230323, "WES2", "/data/project/Meningioma/02.Align/hg38/Dura/05.Final_bam/230323_2_Dura.bam", "chr22", 29661335  )
# Dura_230323_WES1.pysam_read()
# Dura_230323_WES2.pysam_read()
        

# Dura_230405_WES1 = Variant ( 230405, "WES1", "/data/project/Meningioma/02.Align/hg38/Dura/05.Final_bam/230405_Dura.bam", "chr9", 107487067  )
# Dura_230405_WES2 = Variant ( 230405, "WES2", "/data/project/Meningioma/02.Align/hg38/Dura/05.Final_bam/230405_2_Dura.bam", "chr9", 107487067  )
# Dura_230405_duplicated = Variant ( 230405, "duplicated", "/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/01.Pre_bam/230405.sorted.bam", "chr9", 107487067  )
# Dura_230405_deduplicated = Variant ( 230405, "deduplicated", "/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/05.Final_bam/230405_Dura.bam", "chr9", 107487067  )
# Dura_230405_WES1.pysam_read()
# Dura_230405_WES2.pysam_read()
# Dura_230405_duplicated.pysam_read()
# Dura_230405_deduplicated.pysam_read()

# Dura_230419_WES = Variant ( 230419, "WES", "/data/project/Meningioma/02.Align/hg38/Dura/05.Final_bam/230419_Dura.bam", "chr22", 29674837 )
# Dura_230419_duplicated = Variant ( 230419, "duplicated", "/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/01.Pre_bam/230419.sorted.bam", "chr22", 29674837 )
#Dura_230419_deduplicated = Variant ( 230419, "deduplicated", "/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/06.Final_bam/230419_Dura.bam", "chr22", 29674837 )
# Dura_230419_WES.pysam_read()
# Dura_230419_duplicated.pysam_read()
#Dura_230419_deduplicated.pysam_read()
# Dura_230419_duplicated = Variant ( 230419, "duplicated", "/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/01.Pre_bam/230419.sorted.bam", "chr22", 29674836 )
# Dura_230419_duplicated.pysam_read()
# Dura_230419_duplicated = Variant ( 230419, "duplicated", "/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/01.Pre_bam/230419.sorted.bam", "chr22", 29674838 )
# Dura_230419_duplicated.pysam_read()

# Dura_230822_WES = Variant ( 230822, "WES", "/data/project/Meningioma/02.Align/hg38/Dura/05.Final_bam/230822_Dura.bam", "chr22", 29673365  )
# Dura_230822_duplicated = Variant ( 230822, "duplicated", "/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/01.Pre_bam/230822.sorted.bam", "chr22", 29673365  )
# Dura_230822_deduplicated = Variant ( 230822, "deduplicated", "/data/project/Meningioma/00.Targetted/02.Align/hg38/Dura/05.Final_bam/230822_Dura.bam", "chr22", 29673365  )
# Dura_230822_WES.pysam_read()
# Dura_230822_duplicated.pysam_read()
# Dura_230822_deduplicated.pysam_read()