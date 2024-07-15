import argparse
import numpy as np
import pandas as pd
import random
import os, glob, copy, time
import pysam

start_time = time.time()

parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--TISSUE', type=str, default="230822_NF2")
parser.add_argument('--BAMTYPE', type=str, default="Amplicon_single")  #Amplicon_multiplex, Amplicon_single

args = parser.parse_args()

kwargs = {}
kwargs["TISSUE"] = args.TISSUE
kwargs["BAMTYPE"] = args.BAMTYPE


if kwargs["BAMTYPE"] == "Amplicon_multiplex":
    kwargs ["THRESHOLD_END"], kwargs["THRESHOLD_BQ"], kwargs["THRESHOLD_MULTIALLELIC"] = 4, 20, 10000
    kwargs["BAM_DIR"] = "/data/project/Meningioma/51.Amplicon/02.multiplex/02.Align/hg38/01.Pre_bam"
    BAM_PATH_LIST = sorted ( glob.glob ( kwargs["BAM_DIR"]  + "/*.sorted.bam" ) )
    interval_df = pd.read_csv ("/data/project/Meningioma/07.pysam/interval.bed", header = None, sep = "\t")
    output_file = open ("/data/project/Meningioma/51.Amplicon/01.single/07.pysam/log_multiplex.tsv", "w")

elif kwargs["BAMTYPE"] == "Amplicon_single":
    kwargs ["THRESHOLD_END"], kwargs["THRESHOLD_BQ"], kwargs["THRESHOLD_MULTIALLELIC"] = 4, 20, 10000
    kwargs["BAM_DIR"] = "/data/project/Meningioma/51.Amplicon/01.single/02.Align/hg38"

    BAM_PATH_LIST = sorted ( glob.glob ( kwargs["BAM_DIR"] + "/" + kwargs["TISSUE"] + "/01.Pre_bam/*.sorted.bam"  ) ) 
    #BAM_PATH_LIST = ["/data/project/Meningioma/51.Amplicon/01.single/02.Align/hg38/Dura_KLF4/01.Pre_bam/230303_1.sorted.bam", "/data/project/Meningioma/51.Amplicon/01.single/02.Align/hg38/Dura_KLF4/01.Pre_bam/230303_2.sorted.bam", "/data/project/Meningioma/51.Amplicon/01.single/02.Align/hg38/Dura_KLF4/01.Pre_bam/230303_3.sorted.bam" ]
    interval_df = pd.read_csv ("/data/project/Meningioma/07.pysam/interval_{}.bed".format( kwargs["TISSUE"]  ) , header = None, sep = "\t")
    output_file = open ("/data/project/Meningioma/51.Amplicon/02.multiplex/07.pysam/log_{}.tsv".format(  kwargs["TISSUE"]  ), "w")




class Variant:
    def __init__ (self, date, bamtype, BAM_Dir, chr_, pos_, depth_, alt_):
        self.date = date
        self.bamtype = bamtype
        self.BAM_Dir =  BAM_Dir
        self.chr_ = chr_
        self.pos_ = pos_
        self.base_counts = []
        self.depth_ = depth_
        self.alt_ = alt_
        self.reference_base = pysam.Fastafile( "/home/goldpm1/reference/genome.fa" ).fetch(chr_, pos_ - 1, pos_)

   
    def pysam_read (self, **kwargs):
        samfile = pysam.AlignmentFile(self.BAM_Dir, 'rb')
        #for ID in _variant_dic:
        multiallele_check = set()
        ins_count, ins_list, del_count, del_list = 0, [], 0, []
        base_counts = {"A": 0, "T": 0, "C": 0, "G": 0, "N" : 0, "Ins" : 0 , "Del": 0}
        strand_per_base = {"A" : {"F" : 0, "R" : 0},  "T" : {"F" : 0, "R" : 0}, "C" : {"F" : 0, "R" : 0}, "G" : {"F" : 0, "R" : 0}, "N" : {"F" : 0, "R" : 0}, "Ins" : {"F" : 0, "R" : 0}, "Del" : {"F" : 0, "R" : 0}}
        
        iter = samfile.fetch(contig = self.chr_, start = self.pos_-1, end = self.pos_)

        for i, read in enumerate(iter):
            if read.is_proper_pair:
                if read.flag > 512:
                    continue
                read_info = str(read).rstrip().split('\t')
                cigar_string = read.cigarstring
                read_start_pos = read.reference_start + 1

                #1. Soft, hard clip 제거
                # if  ( ('H' not in cigar_string) & ('S' not in cigar_string) ):       
                #     continue

                #2. Read 양쪽 끝에 있으면 제거
                if (self.pos_ >= read.reference_end - kwargs["THRESHOLD_END"] ) | (self.pos_ <= read.reference_start + kwargs ["THRESHOLD_END"] ):  
                    continue
            
                if 'D' in cigar_string:
                    if read_start_pos > self.pos_ :
                        print (cigar_string, self.pos_, read_start_pos)
                    pp = read_start_pos
                    for j in cigar_string.split("D")[0].split("M")[0:-1]:
                        if ("I" in j) | ("S" in j) | ("H" in j):          
                            continue
                        else:
                            pp += int (j)
                    if pp - 1== int(self.pos_):
                        #5. Strand bias 축적
                        if read.is_reverse == False:
                            strand_per_base [ "Del" ]["F"] += 1
                        elif read.is_reverse == True:
                            strand_per_base [ "Del" ]["R"] += 1
                        if  ( ('H' not in cigar_string) & ('S' not in cigar_string) ):       # Soft, hard clip이 아니고 양쪽 끝에 있지 않아야 한다
                            if (self.pos_ < read.reference_end - kwargs["THRESHOLD_END"] ) & (self.pos_ > read.reference_start + kwargs ["THRESHOLD_END"] ):  
                                base_counts["Del"] += 1
                                del_list.append (read.cigarstring)  # read.query_name

                elif 'I' in cigar_string:
                    if read_start_pos > self.pos_ :
                        print (cigar_string, self.pos_, read_start_pos)
                    pp = read_start_pos
                    for j in cigar_string.split("I")[0].split("M")[0:-1]:
                        if ("D" in j) | ("S" in j) | ("H" in j):
                            continue
                        else:
                            pp += int (j)
                    if pp - 1== int(self.pos_):
                        #5. Strand bias 축적
                        if read.is_reverse == False:
                            strand_per_base [ "Ins" ]["F"] += 1
                        elif read.is_reverse == True:
                            strand_per_base [ "Ins" ]["R"] += 1
                        if  ( ('H' not in cigar_string) & ('S' not in cigar_string) ): # Soft, hard clip이 아니고 양쪽 끝에 있지 않아야 한다
                            if (self.pos_ < read.reference_end - kwargs["THRESHOLD_END"] ) & (self.pos_ > read.reference_start + kwargs ["THRESHOLD_END"] ):  
                                base_counts["Ins"] += 1                    
                                ins_list.append (read.cigarstring)

                else:
                    check = {"Clip" : 0,  "End" : 0, "BQ" : 0, "Clustered_event" : 0, "SB" : 0, "Multiallelic" : 0}

                    k = self.pos_ - read_start_pos
                    read_base = read.query_alignment_sequence [ k ]   
                    reference_i = np.where ( np.array(read.get_reference_positions()) == self.pos_ ) [0] [0] - 1
                    reference_base = read.get_reference_sequence() [ reference_i ]                        # 안 맞으면 reference도 소문자로 리턴한다

                    #1. Hard, soft clip 검정
                    if  ( ('H' in cigar_string) | ('S'  in cigar_string) ):       
                        check ["Clip"] = 1

                    #2. End position에 걸려있는지 확인
                    if (self.pos_ <= read.reference_start + kwargs ["THRESHOLD_END"] ) | (self.pos_ >= read.reference_end - kwargs["THRESHOLD_END"] ) :  
                        check ["End"] = 1

                    #3. BQ 검정
                    base_quality = read.query_qualities[ k ]
                    if base_quality < kwargs["THRESHOLD_BQ"]:   #12로 하면 230920을 살리고, 15로 하면 230419를 좋게 만들고
                        check ["BQ"] = 1
                    
                    #4. Clustered event 검정
                    for k2 in range ( max (0, k - 3), k + 4, 1):
                        current_pos = read_start_pos + k2
                        reference_i2 = np.where ( np.array(read.get_reference_positions()) == current_pos ) [0] [0] - 1
                        reference_base2 = read.get_reference_sequence() [ reference_i2 ]
                        read_base2 = read.query_alignment_sequence [ k2 ]
                        if reference_base2 in ["a", "t", "g", "c"]:  # read와 맞지 않는 경우
                            check ["Clustered_event"] += 1
                            
                    #5. Strand bias 축적
                    if read.is_reverse == False:
                        strand_per_base [ read_base.upper() ]["F"] += 1
                    elif read.is_reverse == True:
                        strand_per_base [ read_base.upper() ]["R"] += 1
                        

                    if ( check["Clip"] == 0 ) & (check ["End"] == 0) & (check["BQ"] == 0) & ( check ["Clustered_event"]  in [0, 1] ):   #통과해야만
                        base_counts [ read_base ] += 1
                        # if reference_base != read_base:
                        #     print ( "\ti = {}\tk = {}\tcurrent_pos = {}\treference_base = {}\tread_base = {}\tbase_quality = {}".format (i, k, self.pos_, reference_base, read_base, base_quality ))


        
        print ("\t# of total read depth : {}". format ( np.sum ( list ( base_counts.values() )  ) ) )
        #print ( "\treference_base = {}".format ( self.reference_base ) )
        for key in base_counts.keys():
            print ("\t# of {} : {}\tstrand = {}". format (key, base_counts[key], strand_per_base[key] ))
            if kwargs["BAMTYPE"] != "Amplicon_single":
                if key != self.reference_base:
                    if (strand_per_base[key]["F"] == 0) | (strand_per_base[key]["R"] == 0):    # Strand bias에 걸린 경우
                        if ( strand_per_base[key]["F"] + strand_per_base[key]["R"] >= 3) :  
                            base_counts[key] = 0
                            print ("\t(수정)# of {} : {}\tstrand = {}". format (key, base_counts[key], strand_per_base[key] ))
        #print ("\t\tdel_list = {}".format (del_list))

        # Multiallelic count 위해
        if kwargs["BAMTYPE"] != "Amplicon_single":
            if kwargs["BAMTYPE"] == "Amplicon_multiplex":
                kwargs["THRESHOLD_MULTIALLELIC"] = sum( base_counts.values() ) * 0.02   # 0.02%가 noise의 한계 point라고 생각

            c = 0
            for key in base_counts.keys():
                if key in ["Ins", "Del"]:
                    continue

                if (strand_per_base[key]["F"] >= kwargs["THRESHOLD_MULTIALLELIC"] ) & (strand_per_base[key]["R"] < kwargs["THRESHOLD_MULTIALLELIC"]) & (base_counts[key] >= kwargs["THRESHOLD_MULTIALLELIC"]):
                    c += 1
                    print ( "\t{}\t{}\tthreshold = {}".format (key, strand_per_base[key], kwargs["THRESHOLD_MULTIALLELIC"]))
                elif (strand_per_base[key]["R"] >= kwargs["THRESHOLD_MULTIALLELIC"] ) & (strand_per_base[key]["F"] < kwargs["THRESHOLD_MULTIALLELIC"]) & (base_counts[key] >= kwargs["THRESHOLD_MULTIALLELIC"]) :
                    c += 1
                    print ( "\t{}\t{}\tthreshold = {}".format (key, strand_per_base[key], kwargs["THRESHOLD_MULTIALLELIC"]))
            if c >= 2:
                for key in base_counts.keys():
                    if key != self.reference_base:   # Alt는 죄다 0으로 만들어줌
                        base_counts[key] = 0    
                    print ("\t(수정) # of {} : {}\tstrand = {}". format (key, base_counts[key], strand_per_base[key] ))

        self.base_counts = base_counts
        self.depth_ = sum( base_counts.values() )


def HEATMAP_VISUALIZATION (vaf_df, title, Output_filename, **kwargs):
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcl
    import numpy as np
    import pandas as pd
    from matplotlib.colors import LinearSegmentedColormap

    plt.rcParams["font.family"] = 'arial'

    vaf_df = vaf_df  * 100
    vaf_df = vaf_df.applymap(lambda x: np.floor(x * 10) / 10)
    #vaf_df = vaf_df.round(1)

    # Define the colors
    colors = [ "white", "#A3B18A", "#588157"]
    fig, ax = plt.subplots ( nrows = 1, ncols = 1, figsize =(9 / 2.54, 8 / 2.54))

    positions = [0, 0.01, 1]  # Define the positions for each color
    # Create the colormap
    cmap = LinearSegmentedColormap.from_list('custom_cmap', list(zip(positions, colors)))

    fig.subplots_adjust ( wspace = 0.4, bottom = 0.03, top = 0.7, left = 0.22, right = 0.98)
    fig.set_facecolor('white')

    sns.heatmap (vaf_df , cmap = cmap, linewidths = 0.5, linecolor = "black", annot = vaf_df, annot_kws={"size": 25 / np.sqrt(len(vaf_df))} )   # fmt=".2f", 
    fig.suptitle ( title, fontsize = 12, fontweight = "bold", ha = "left", x = 0 )
    ax.set_xticklabels( ax.get_xticklabels(), fontsize = 7, ha = 'left' )
    ax.tick_params(axis = 'x',  rotation = 45, pad = -2.5)
    ax.set_yticklabels( [ i.get_text().replace( "_Multiplex.sorted" , "" ) for i in ax.get_yticklabels()], fontsize = 5, va = 'center' )
    ax.tick_params(axis = 'y', pad = -1.5 )

    plt.tick_params(axis='both', which='major', labelsize = 7, left = False, labelbottom = False, bottom=False, top = False, labeltop=True)

    fig.savefig ( Output_filename, dpi = 300)


def BAMSNAP_IGV ( BAM_DIR_LIST, TITLE_LIST, POS, OUT) :
    command = " ".join (  ["bamsnap -bam ", BAM_DIR_LIST, "-title", TITLE_LIST, "-pos", POS, "-out", OUT,
                        "-draw coordinates bamplot base gene", 
                        "-bamplot coverage base read",
                        "-refversion hg38", 
                        "-read_group strand",
                        "-margin 10",
                        "-height 50",
                        "-plot_margin_left 20",
                        "-plot_margin_right 20",
                        "-border",
                        "-base_height 50" ]
                        )
    print (command)
    os.system (command)




vaf_df = pd.DataFrame ( np.zeros ( ( len(BAM_PATH_LIST), interval_df.shape[0] ) ))
interval_df ["ID"]  = interval_df.iloc[:, 0].astype(str) + ":" + interval_df.iloc[:, 1].astype(str) + "(" + interval_df.iloc[:, 4].astype(str) + ")"
vaf_df.columns = list (interval_df ["ID"])
count_df = pd.DataFrame ( [['' for _ in range(interval_df.shape[0])] for _ in range( len(BAM_PATH_LIST) )] ) 


row_list = []
count_matrix = []

    

for BAM_PATH_i, BAM_PATH in enumerate ( BAM_PATH_LIST ) :
    DATE = BAM_PATH.split("/")[-1].split(".sorted.bam")[0]

    print ( "\n\n{}".format( DATE )  )
    print (DATE, file = output_file)
    
    row_list.append ( DATE )
    count_list = []

    for k in range (interval_df.shape[0]): 
    #for k in range (7, 8):  # chr9    107487067   T   G   KLF4    SNP 230303,230405
        CHR, POS, REF, ALT, GENE, TYPE = interval_df.iloc[k][0], interval_df.iloc[k][1], interval_df.iloc[k][2], interval_df.iloc[k][3], interval_df.iloc[k][4], interval_df.iloc[k][5]
        ALT_TYPE = TYPE if TYPE in ["Ins", "Del"] else ALT

        print (CHR, POS)
        Count = Variant ( DATE, kwargs["BAMTYPE"], BAM_PATH, CHR, POS, 0, ALT_TYPE  )
        print ( "\t{}:{}\t{}>{}".format (Count.chr_, Count.pos_, REF, ALT), end = "\t")
        Count.pysam_read ( **kwargs )

        vaf_df.iloc [BAM_PATH_i, k] =  Count.base_counts [Count.alt_]  /  Count.depth_  if Count.depth_ != 0 else 0
        count_list.append ( "{}/{}".format( str(Count.base_counts [Count.alt_])  ,  str(Count.depth_)  ) )
        print ( "\t{}:{} {}>{}\t{}% ( = {} / {})\t{}".format (Count.chr_, Count.pos_, REF, ALT, round(vaf_df.iloc [BAM_PATH_i, k] * 100, 2), Count.base_counts [Count.alt_] , Count.depth_, Count.base_counts) )
        print ( "\t{}:{} {}>{}\t{}% ( = {} / {})\t{}".format (Count.chr_, Count.pos_, REF, ALT, round(vaf_df.iloc [BAM_PATH_i, k] * 100, 2), Count.base_counts [Count.alt_] , Count.depth_, Count.base_counts) , file = output_file)
    print("Time elapsed : {}min".format ( round ( (time.time() - start_time) / 60 ) ) )

    count_matrix.append (count_list)
output_file.close()

vaf_df.index = row_list

count_df = pd.DataFrame (count_matrix, columns = vaf_df.columns )
count_df.index = vaf_df.index


if kwargs ["BAMTYPE"] == "Amplicon_single":
    vaf_df.to_csv ("/data/project/Meningioma/51.Amplicon/01.single/07.pysam/{}_vaf_df.tsv".format(kwargs["TISSUE"]), sep = "\t")
    count_df.to_csv ("/data/project/Meningioma/51.Amplicon/01.single/07.pysam/{}_count_df.tsv".format(kwargs["TISSUE"]), sep = "\t")
    HEATMAP_VISUALIZATION (vaf_df, "Amplicon_".format (kwargs["TISSUE"]), "/data/project/Meningioma/51.Amplicon/01.single/07.pysam/{}_heatmap.pdf".format (kwargs["TISSUE"]), **kwargs )
elif kwargs ["BAMTYPE"] == "Amplicon_multiplex":
    vaf_df.to_csv ("/data/project/Meningioma/51.Amplicon/02.multiplex/07.pysam/multiplex_vaf_df.tsv", sep = "\t")
    count_df.to_csv ("/data/project/Meningioma/51.Amplicon/02.multiplex/07.pysam/multiplex_count_df.tsv", sep = "\t")
    HEATMAP_VISUALIZATION (vaf_df, "Amplicon_multiplex", "/data/project/Meningioma/51.Amplicon/02.multiplex/07.pysam/multiplex_heatmap.pdf", **kwargs )


