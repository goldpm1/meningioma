import os, subprocess, argparse
import pandas as pd
import numpy as np
import vcf, pysam
from tqdm import tqdm
from scipy.stats import binom, betabinom

kwargs = {}

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument('--INPUT_VCF', type=str)
parser.add_argument('--OUTPUT_VCF', type=str)
parser.add_argument('--THRESHOLD_DEPTH', type=int)
parser.add_argument('--OUTPUT_BAMSNAP_DIR', type=str)
parser.add_argument('--BAM_DIR_LIST', type=str)
parser.add_argument('--TITLE_LIST', type=str)
parser.add_argument('--SAMPLE_ID', type=str)
parser.add_argument('--WES_TUMOR_BED', type=str)
parser.add_argument('--WES_DURA_VCF', type=str)

args = parser.parse_args()

kwargs["INPUT_VCF"] = args.INPUT_VCF
kwargs["OUTPUT_VCF"] = args.OUTPUT_VCF
kwargs["THRESHOLD_DEPTH"] = args.THRESHOLD_DEPTH
kwargs["OUTPUT_BAMSNAP_DIR"] = args.OUTPUT_BAMSNAP_DIR
kwargs["BAM_DIR_LIST"] = args.BAM_DIR_LIST
kwargs["TITLE_LIST"] = args.TITLE_LIST
kwargs["SAMPLE_ID"] = args.SAMPLE_ID
kwargs["WES_TUMOR_BED"] = args.WES_TUMOR_BED
kwargs["WES_DURA_VCF"] = args.WES_DURA_VCF



for subdir in ["01.Image_singlecell", "02.Image_WES_tumor",  "03.Image_WES_dura"]:
    if not os.path.exists ( kwargs["OUTPUT_BAMSNAP_DIR"] + "/" + subdir ):
        os.makedirs (  kwargs["OUTPUT_BAMSNAP_DIR"] + "/" + subdir )



####################################################################################################################################

kwargs["BAMTYPE"] = "WES"
kwargs ["THRESHOLD_END"], kwargs["THRESHOLD_BQ"], kwargs["THRESHOLD_MULTIALLELIC"] = 4, 15, 2


class Variant:
    def __init__ (self, date, bamtype, BAM_Dir, chr_, pos_):
        self.date = date
        self.bamtype = bamtype
        self.BAM_Dir =  BAM_Dir
        self.chr_ = chr_
        self.pos_ = pos_
        self.base_counts = []
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
                    for k2 in range ( max (0, k - (kwargs ["THRESHOLD_END"] - 1 )), k + ( kwargs ["THRESHOLD_END"] ), 1):
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
                        #     print ( "\tk = {}\tcurrent_pos = {}\treference_base = {}\tread_base = {}\tbase_quality = {}".format (k, self.pos_, reference_base, read_base, base_quality ))


        
        # print ("\t# of total read depth : {}". format ( np.sum ( list ( base_counts.values() )  ) ) )
        # print ( "\treference_base = {}".format ( self.reference_base ) )
        for key in base_counts.keys():
            #print ("\t# of {} : {}\tstrand = {}". format (key, base_counts[key], strand_per_base[key] ))
            if key != self.reference_base:
                if (strand_per_base[key]["F"] == 0) | (strand_per_base[key]["R"] == 0):    # Strand bias에 걸린 경우
                    if ( strand_per_base[key]["F"] + strand_per_base[key]["R"] >= 3) :  
                        base_counts[key] = 0
                        #print ("\t(수정)# of {} : {}\tstrand = {}". format (key, base_counts[key], strand_per_base[key] ))
        #print ("\t\tdel_list = {}".format (del_list))

        # Multiallelic count 위해
        c = 0
        for key in base_counts.keys():
            if (strand_per_base[key]["F"] >= kwargs["THRESHOLD_MULTIALLELIC"] ) & (strand_per_base[key]["R"] < kwargs["THRESHOLD_MULTIALLELIC"]) & (base_counts[key] >= kwargs["THRESHOLD_MULTIALLELIC"]):
                c += 1
                #print ( "\t{}\t{}\t{}".format (key, strand_per_base[key], kwargs["THRESHOLD_MULTIALLELIC"]))
            elif (strand_per_base[key]["R"] >= kwargs["THRESHOLD_MULTIALLELIC"] ) & (strand_per_base[key]["F"] < kwargs["THRESHOLD_MULTIALLELIC"]) & (base_counts[key] >= kwargs["THRESHOLD_MULTIALLELIC"]) :
                c += 1
                #print ( "\t{}\t{}\t{}".format (key, strand_per_base[key], kwargs["THRESHOLD_MULTIALLELIC"]))
        if c >= 2:
            for key in base_counts.keys():
                if key != self.reference_base:   # Alt는 죄다 0으로 만들어줌
                    base_counts[key] = 0    
                #print ("\t(수정) # of {} : {}\tstrand = {}". format (key, base_counts[key], strand_per_base[key] ))

        self.base_counts = base_counts
        self.depth_ = sum( base_counts.values() )




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
    #print (command)
    os.system (command)


def print_line_vcf (line, rescue_list):
    # Print vcf
    print_line = [ str( line.CHROM) , str( line.POS ) , ".", str (line.REF ), str (line.ALT[0] ), str( line.QUAL ) , "." ]
    temp = []
    for key, value in line.INFO.items():
        temp.append ( f"{key}={value}" )
    print_line.append ( ';'.join (temp).replace(" ", "").replace("\'", "").replace("[", "").replace("]", "") )
    print_line.append ( "GT:AD:DP:GQ:PL" )

    # Print information for each sample
    for sample in line.samples:
        GT = str(sample.data.GT).replace(" ", "").replace("[", "").replace("]", "")
        if sample.sample in rescue_list:
            GT = "0/1"
        

        print_line.append ( ":".join (  [ GT, 
                                                            str(sample.data.AD).replace(" ", "").replace("[", "").replace("]", ""), 
                                                            str(sample.data.DP).replace(" ", "").replace("[", "").replace("]", "").replace("None", "."),
                                                            str(sample.data.GQ).replace(" ", "").replace("[", "").replace("]", "").replace("None", "."), 
                                                            str(sample.data.PL).replace(" ", "").replace("[", "").replace("]", "").replace("None", ".")  ] 
                                                            ) )

    return print_line

def find_blood_sample(samplenames):
  for i, name in enumerate(samplenames):
    if "Blood" in name:
      return i
  return -1

def binom_filter ( depth, alt , p ):
    return binom.cdf( alt, depth, p)

def betabinom_filter ( depth, alt , a, b ):
    return betabinom.cdf( alt, depth, a, b )



#####################################################################################################################################


# 읽고쓰기
input_file = open(kwargs["INPUT_VCF"], "r")
output_file = open ( kwargs ["OUTPUT_VCF"], "w")
for line in input_file.readlines():
    if line[0] != '#':
        break
    print (line.rstrip(), file = output_file)
input_file.close()
output_file.close()

# os.system ("rm -rf " + kwargs["OUTPUT_BAMSNAP_DIR"] +  "/01.Image_singlecell/*")
# os.system ("rm -rf " + kwargs["OUTPUT_BAMSNAP_DIR"] +  "/02.Image_WES_tumor/*")
# os.system ("rm -rf " + kwargs["OUTPUT_BAMSNAP_DIR"] +  "/03.Image_WES_dura/*")

input_file = open ( kwargs["INPUT_VCF"], "r")
output_file = open ( kwargs ["OUTPUT_VCF"], "w")
# Print header
for line in input_file:
    line = line.rstrip()
    if line[0] != "#":
        break
    print (line, file = output_file)
input_file.close()

vcf_reader = vcf.Reader(open( kwargs["INPUT_VCF"], "r"))


cnt = 0
#for line in tqdm ( vcf_reader ) :        # line.CHROM, recrod.POS ,line.ALT
for line in vcf_reader :        # line.CHROM, recrod.POS ,line.ALT
    if ( len(line.ALT) >= 2 ):  # multiallelic 은 빼자
        continue

    CHR = line.CHROM
    POS = str (line.POS)
    REF, ALT = str(line.REF), str( line.ALT[0] )
    samplenames = vcf_reader.samples
    sample_i_Blood = find_blood_sample ( samplenames )

    # if (CHR != "chr2") | (line.POS < 237372266):
    #     continue

    # Blood를 처리
    if ( line.samples [ sample_i_Blood ].data.DP  == None ):
        continue
    if int( line.samples [ sample_i_Blood ].data.DP ) < kwargs["THRESHOLD_DEPTH"]:
        continue
    if ( line.samples [sample_i_Blood].data.GT in ["0/0", "0|0"] ):                # Blood가 진짜 0/0, 0|0 인지 binomial filter를 통해 확인 (진짜 alt, depth를 센다)
        ALT_TYPE = "Ins" if len(REF) < len (ALT) else "Del" if len(REF) > len(ALT) else ALT
        Count_blood = Variant ( kwargs["TITLE_LIST"].split(",") [ -1 ], "WES", kwargs["BAM_DIR_LIST"].split(",") [ -1 ], CHR, int(POS)   )
        Count_blood.pysam_read ( **kwargs )
        ALT_COUNTS, DEPTH_COUNTS = Count_blood.base_counts [ ALT_TYPE ],  Count_blood.depth_
        binom_p_blood = binom_filter (  DEPTH_COUNTS,  ALT_COUNTS, 0.05  )
        if binom_p_blood > 0.97:
            # print ("\n{}:{}\tALT_COUNTS = {}\tDEPTH_COUNTS = {}\tp = {}\t(HC)BloodGT: {} -> 0/1로 변경해야 할듯".format (  line.CHROM, line.POS,  ALT_COUNTS, DEPTH_COUNTS, binom_p_blood,line.samples [sample_i_Blood].data.GT  ))
            # print ( Count_blood.base_counts  )
            check = 1
            continue


    # 다른 sample들을 대상으로
    check = 0
    ALT_COUNTS_LIST = np.zeros ( len (samplenames) )
    DEPTH_COUNTS_LIST = np.zeros ( len (samplenames) )

    for sample_i, samplename in enumerate (samplenames):
        confirm_list = []
        if sample_i == sample_i_Blood:         # Blood sample이 아닌 것을 대상으로
            continue

        # 일단 depth, alt를 무조건 다시 계산
        ALT_TYPE = "Ins" if len(REF) < len (ALT) else "Del" if len(REF) > len(ALT) else ALT
        Count = Variant ( kwargs["TITLE_LIST"].split(",") [ sample_i ], "WES", kwargs["BAM_DIR_LIST"].split(",") [ sample_i ], CHR, int(POS)  )
        Count.pysam_read ( **kwargs )
        ALT_COUNTS, DEPTH_COUNTS = Count.base_counts [ ALT_TYPE ],  Count.depth_
        ALT_COUNTS_LIST [sample_i], DEPTH_COUNTS_LIST [sample_i] = int(ALT_COUNTS), int(DEPTH_COUNTS)


        if line.samples [sample_i].data.DP  == None:
            continue
        if int( line.samples [sample_i].data.DP ) < kwargs["THRESHOLD_DEPTH"]:   # Depth가 너무 낮으면 pass
            continue

        #print ("\n{}\t{}:{}\tdepth = {}\talt = {}".format (samplename, line.CHROM, line.POS, DEPTH_COUNTS, ALT_COUNTS ))

        if  ( line.samples [sample_i].data.GT in ["1/1", "1|1"] ) & ( line.samples [sample_i_Blood].data.GT in ["0/0", "0|0"] ):  # Homo이고 Blood는 아닐 경우
            binom_p_sample = binom_filter (  DEPTH_COUNTS,  ALT_COUNTS, 0.9  )
            binom_p_blood = binom_filter (  int (line.samples [sample_i_Blood].data.DP),  int (line.samples [sample_i_Blood].data.AD[1]), 0.01  )           # 0.01은 너무 빡셀 것 같은데...

            if (binom_p_sample > 0.05) :   # 진짜 1|1인지 binomial filter
                if (binom_p_blood < 0.97):  # 진짜 0|0 인지 binomial filter
                    check = 1; confirm_list.append ( samplename )
                    print ("{}\t{}:{}\tdepth = {}\talt = {}\tSampleGT: {} (binom p = {})\tBloodGT: {} (p = {})".format (samplename, line.CHROM, line.POS, DEPTH_COUNTS, ALT_COUNTS,
                                                                            line.samples [sample_i].data.GT, round (binom_p_sample, 2) ,  line.samples [sample_i_Blood].data.GT, round (binom_p_blood, 2)  ))

        if  ( line.samples [sample_i].data.GT in ["0/1", "0|1"] ) & ( line.samples [sample_i_Blood].data.GT in ["0/0", "0|0"] ):  # Hetero이고 Blood는 아닐 경우
            binom_p_sample = binom_filter (  DEPTH_COUNTS,  ALT_COUNTS, 0.5  )
            binom_p_blood = binom_filter (  int (line.samples [sample_i_Blood].data.DP),  int (line.samples [sample_i_Blood].data.AD[1]), 0.01  )

            # print ("\n{}\t{}:{}\tdepth = {}\talt = {}\tSampleGT: {} (p = {})\tBloodGT: {} (p = {})".format (samplename, line.CHROM, line.POS, DEPTH_COUNTS, ALT_COUNTS,
            #                                                                                                 line.samples [sample_i].data.GT, round (binom_p_sample, 2) ,  line.samples [sample_i_Blood].data.GT, round (binom_p_blood, 2)  ))

            if (binom_p_sample > 0.05) & (binom_p_sample < 0.97):   # 진짜 0|1인지 binomial filter
                if (binom_p_blood < 0.97):  # 진짜 0|0 인지 binomial filter
                    check = 3;  confirm_list.append ( samplename )
                    print ("{}\t{}:{}\tdepth = {}\talt = {}\tSampleGT: {} (binom p = {})\tBloodGT: {} (p = {})".format (samplename, line.CHROM, line.POS, DEPTH_COUNTS, ALT_COUNTS,
                                                                                                line.samples [sample_i].data.GT, round (binom_p_sample, 2) ,  line.samples [sample_i_Blood].data.GT, round (binom_p_blood, 2)  ))
                    



    # 더 완화된 betabinomial을 기준으로 해서 다시한번 기회를 줌
    if check != 0:
        rescue_list, rescue_depth_list = [], []
        for sample_i, samplename in enumerate (samplenames):
            ALT_COUNTS, DEPTH_COUNTS = int ( ALT_COUNTS_LIST [sample_i] ) , int ( DEPTH_COUNTS_LIST [sample_i] )

            if sample_i == sample_i_Blood:         # Blood sample이 아닌 것을 대상으로
                continue
            if ( line.samples [sample_i].data.DP  == None )  :
                check = 99; rescue_list.append ( samplename); rescue_depth_list.append ( samplename )
                print ("\t→ {}\tRESCUE 확정 (그냥 살려줌): depth = {}\talt = {}\tSampleGT: {} → 0/1".format (samplename, DEPTH_COUNTS, ALT_COUNTS,
                                                                                                line.samples [sample_i].data.GT  )  )
                continue
            elif  (int( line.samples [sample_i].data.DP ) < kwargs["THRESHOLD_DEPTH"] / 3) :
                check = 99; rescue_list.append ( samplename); rescue_depth_list.append ( samplename )
                print ("\t→ {}\tRESCUE 확정 (그냥 살려줌): depth = {}\talt = {}\tSampleGT: {} → 0/1".format (samplename, DEPTH_COUNTS, ALT_COUNTS,
                                                                                                line.samples [sample_i].data.GT  )  )
                continue

        
            binom_p_sample = binom_filter (  DEPTH_COUNTS,  ALT_COUNTS, 0.5  )
            betabinom_p_sample = betabinom_filter (  DEPTH_COUNTS,  ALT_COUNTS, DEPTH_COUNTS / 5,  DEPTH_COUNTS / 5  )
            binom_p_blood = binom_filter (  int (line.samples [sample_i_Blood].data.DP),  int (line.samples [sample_i_Blood].data.AD[1]), 0.01  )

            # 아까 그놈이면 가라
            if (binom_p_sample > 0.05) & (binom_p_sample < 0.97):   # 진짜 0|1인지 binomial filter
                if (binom_p_blood < 0.97):  # 진짜 0|0 인지 binomial filter
                    check = 3
                    continue
            if  ( line.samples [sample_i].data.GT in ["0/1", "0|1", "1/1", "1|1"] ):
                continue



            # Resuce 심사대상
            print ("\t{}\t{}:{}\tRECUE 심사대상: depth = {}\talt = {}\tSampleGT: {} (betabinom p = {})\tBloodGT: {} (binom p = {})".format (samplename, line.CHROM, line.POS, DEPTH_COUNTS, ALT_COUNTS,
                                                                                                            line.samples [sample_i].data.GT, round (binom_p_sample, 2) ,  line.samples [sample_i_Blood].data.GT, round (binom_p_blood, 2)  ))

            if ( betabinom_p_sample >= 0.02) & (binom_p_sample < 0.97):   # 진짜 0|1인지 binomial filter
                if (binom_p_blood < 0.97):  # 진짜 0|0 인지 binomial filter
                    check = 99; rescue_list.append ( samplename)   
                    print ("\t→ {}\tRESCUE 확정: depth = {}\talt = {}\tSampleGT: {} (betabinom p = {})\tBloodGT: {} (binom p = {})".format (samplename, DEPTH_COUNTS, ALT_COUNTS,
                                                                                                line.samples [sample_i].data.GT, round (betabinom_p_sample, 2) ,  line.samples [sample_i_Blood].data.GT, round (binom_p_blood, 2)  ))
            # else:
            #     print ("\t→ {}\tRESCUE 탈락: depth = {}\talt = {}\tSampleGT: {} (betabinom p = {})".format (samplename, DEPTH_COUNTS, ALT_COUNTS,
            #                                                                                     line.samples [sample_i].data.GT, round (betabinom_p_sample, 2)   ))
        
                
        if len (rescue_list) > len (samplename) / 2:
            print ("\t너무 많이 살려줘야 해서 그냥 다 기각")
            rescue_list = []
        if len (rescue_depth_list) > len (confirm_list):
            print ("\t너무 많이 살려줘야 해서 (>confirm_list) 그냥 다 기각")
            rescue_list = []

        # 한줄씩 출력
        print ("\trescue list : {}".format (rescue_list))
        print_line = print_line_vcf ( line, rescue_list )
        print ("\t".join ( print_line), file = output_file )
        
        cnt += 1
        print ("")        

        # ########### BAMSNAP도 찍자 (conda activate cnvpytor) ######################
        # BAMSNAP_IGV ( kwargs ["BAM_DIR_LIST"].replace(",", " "),
        #                             kwargs ["TITLE_LIST"].replace(",", " "), 
        #                             "{}:{}".format (CHR, str(POS)) , 
        #                             kwargs["OUTPUT_BAMSNAP_DIR"] +  "/01.Image_singlecell/{}_{}:{}.jpg".format( kwargs["SAMPLE_ID"] , CHR, str(POS) )) 

output_file.close()

os.system ( "bgzip -c -f {} > {}.gz".format ( kwargs ["OUTPUT_VCF"], kwargs ["OUTPUT_VCF"] ))
os.system ( "tabix -p vcf {}.gz".format (kwargs ["OUTPUT_VCF"]) )


#### WES Tumor에서 나왔던 position도 그림으로 찍어보기
df_interval = pd.read_csv( kwargs["WES_TUMOR_BED"], sep = "\t" )
for k in range ( df_interval.shape[0] ):
    CHR = df_interval.iloc[k][0]
    POS = str( df_interval.iloc[k][2] )

    BAMSNAP_IGV ( kwargs ["BAM_DIR_LIST"].replace(",", " "),
                            kwargs ["TITLE_LIST"].replace(",", " "), 
                            "{}:{}".format (CHR, str(POS)) , 
                            kwargs["OUTPUT_BAMSNAP_DIR"] +  "/02.Image_WES_tumor/{}_{}:{}.jpg".format( kwargs["SAMPLE_ID"] , CHR, str(POS) )) 


#### WES Durar에서 나왔던 position도 그림으로 찍어보기
vcf_dura = vcf.Reader(  open( kwargs["WES_DURA_VCF"], "r") )
for line in vcf_dura:        # line.CHROM, recrod.POS ,line.ALT
    CHR = line.CHROM
    POS = str (line.POS)
    
    BAMSNAP_IGV ( kwargs ["BAM_DIR_LIST"].replace(",", " "),
                            kwargs ["TITLE_LIST"].replace(",", " "), 
                            "{}:{}".format (CHR, str(POS)) , 
                            kwargs["OUTPUT_BAMSNAP_DIR"] +  "/03.Image_WES_dura/{}_{}:{}.jpg".format( kwargs["SAMPLE_ID"] , CHR, str(POS) )) 
