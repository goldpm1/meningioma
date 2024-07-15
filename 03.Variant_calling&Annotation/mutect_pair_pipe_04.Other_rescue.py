import argparse
import numpy as np
import pandas as pd
import os, pysam, vcf, pybedtools
import gzip

# python3 /data/project/Meningioma/script/03.Variant_calling&Annotation/mutect_pair_pipe_04.Other_rescue.py

parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--Sample_ID', type=str, default = "230419")
parser.add_argument('--TISSUE', type=str, default = "Dura")
parser.add_argument('--MINIMUM_ALT', type=int, default = 1)
parser.add_argument('--TUMOR_INTERVAL', type=str, default = "/data/project/Meningioma/04.mutect/03.Tumor_interval/230419_Tumor.MT2.FMC.HF.RMBLACK.bed")
parser.add_argument('--CASE_BAM_PATH', type=str, default = "/data/project/Meningioma/02.Align/hg38/Dura/05.Final_bam/230419_Dura.bam")
parser.add_argument('--CONTROL_BAM_PATH', type=str, default = "/data/project/Meningioma/02.Align/hg38/Blood/05.Final_bam/230419_Blood.bam")
parser.add_argument('--TUMOR_MUTECT2_VCF', type=str, default = "/data/project/Meningioma/04.mutect/02.PASS/230419_Tumor.MT2.FMC.HF.RMBLACK.vcf")
parser.add_argument('--OTHER_MUTECT2_VCF', type=str, default = "/data/project/Meningioma/04.mutect/02.PASS/230419_Dura.MT2.FMC.HF.RMBLACK.vcf")
parser.add_argument('--HC_GVCF', type=str, default = "/data/project/Meningioma/05.gvcf/02.remove_nonref/230419/230419_Dura.g.vcf")
parser.add_argument('--RESCUE_VCF', type=str, default = "/data/project/Meningioma/04.mutect/04.Other_rescue/230419_Dura.MT2.FMC.HF.RMBLACK.rescue.vcf")
parser.add_argument('--TUMOR_SHARED_VARIANT_VCF', type=str, default = "/data/project/Meningioma/04.mutect/05.Shared_variant/230419_Tumor.MT2.FMC.HF.RMBLACK.shared_variant.vcf")
parser.add_argument('--OTHER_SHARED_VARIANT_VCF', type=str, default = "/data/project/Meningioma/04.mutect/05.Shared_variant/230419_Dura.MT2.FMC.HF.RMBLACK.shared_variant.vcf")
parser.add_argument('--TUMOR_UNIQUE_VCF', type=str, default = "/data/project/Meningioma/04.mutect/06.Unique/230419_Dura.MT2.FMC.HF.RMBLACK.tumor_unique.vcf")
parser.add_argument('--OTHER_UNIQUE_VCF', type=str, default = "/data/project/Meningioma/04.mutect/06.Unique/230419_Dura.MT2.FMC.HF.RMBLACK.other_unique.vcf")

args = parser.parse_args()

kwargs = {}
kwargs["Sample_ID"] = args.Sample_ID
kwargs["TISSUE"] = args.TISSUE
kwargs["MINIMUM_ALT"] = args.MINIMUM_ALT
kwargs["TUMOR_INTERVAL"] = args.TUMOR_INTERVAL
kwargs["CASE_BAM_PATH"] = args.CASE_BAM_PATH
kwargs["CONTROL_BAM_PATH"] = args.CONTROL_BAM_PATH
kwargs["TUMOR_MUTECT2_VCF"] = args.TUMOR_MUTECT2_VCF
kwargs["OTHER_MUTECT2_VCF"] = args.OTHER_MUTECT2_VCF
kwargs["HC_GVCF"] = args.HC_GVCF
kwargs["OTHER_RESCUE_VCF"] = args.RESCUE_VCF
kwargs["TUMOR_SHARED_VARIANT_VCF"] = args.TUMOR_SHARED_VARIANT_VCF
kwargs["OTHER_SHARED_VARIANT_VCF"] = args.OTHER_SHARED_VARIANT_VCF
kwargs["TUMOR_UNIQUE_VCF"] = args.TUMOR_UNIQUE_VCF
kwargs["OTHER_UNIQUE_VCF"] = args.OTHER_UNIQUE_VCF

kwargs ["THRESHOLD_END"], kwargs["THRESHOLD_BQ"], kwargs["THRESHOLD_MULTIALLELIC"] = 4, 15, 1

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
                        if ("I" in j) | ("S" in j):          
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
                        if ("D" in j) | ("S" in j):
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
                        #     print ( "\tk = {}\tcurrent_pos = {}\treference_base = {}\tread_base = {}\tbase_quality = {}".format (k, self.pos_, reference_base, read_base, base_quality ))


        if self.bamtype == "Dura":
            print ("\t# of total read depth : {}". format ( np.sum ( list ( base_counts.values() )  ) ) )
    
        for key in base_counts.keys():
            if self.bamtype == "Dura":
                print ("\t# of {} : {}\tstrand = {}". format (key, base_counts[key], strand_per_base[key] ))
            if key != self.reference_base:
                if (strand_per_base[key]["F"] == 0) | (strand_per_base[key]["R"] == 0):    # Strand bias에 걸린 경우
                    if ( strand_per_base[key]["F"] + strand_per_base[key]["R"] >= 3) :  
                        base_counts[key] = 0
                        if self.bamtype == "Dura":
                            print ("\t(Strand bias 수정)# of {} : {}\tstrand = {}". format (key, base_counts[key], strand_per_base[key] ))
        if self.bamtype == "Dura":
            print ("\t\tdel_list = {}".format (del_list))

        # Multiallelic count 위해
        c = 0
        for key in base_counts.keys():
            if (strand_per_base[key]["F"] >= kwargs["THRESHOLD_MULTIALLELIC"] ) & (strand_per_base[key]["R"] < kwargs["THRESHOLD_MULTIALLELIC"]):
                c += 1
            elif (strand_per_base[key]["R"] >= kwargs["THRESHOLD_MULTIALLELIC"] ) & (strand_per_base[key]["F"] < kwargs["THRESHOLD_MULTIALLELIC"]):
                c += 1
        if c >= 2:
            for key in base_counts.keys():
                if key != self.reference_base:   # Alt는 죄다 0으로 만들어줌
                    base_counts[key] = 0    
                if self.bamtype == "Dura":
                    print ("\t(Multiallelic수정) # of {} : {}\tstrand = {}". format (key, base_counts[key], strand_per_base[key] ))

        self.base_counts = base_counts
        self.depth_ = sum( base_counts.values() )


def print_line_vcf (line):
    # Print vcf
    print_line = [ line.CHROM, line.POS, ".", line.REF, line.ALT[0], ".", line.FILTER ]
    temp = []
    for key, value in line.INFO.items():
        temp.append ( f"{key}={value}" )
    print_line.append ( ';'.join (temp).replace(" ", "").replace("\'", "").replace("[", "").replace("]", "") )
    print_line.append ( "GT:AD:AF:DP:F1R2:F2R1:FAD:SB" )

    # Print information for each sampl
    for sample in line.samples:
        print_line.append ( ":".join (  [str(sample.data.GT).replace(" ", "").replace("\'", "").replace("[", "").replace("]", ""), 
                                                            str(sample.data.AD).replace(" ", "").replace("\'", "").replace("[", "").replace("]", ""), 
                                                            str(sample.data.AF).replace(" ", "").replace("\'", "").replace("[", "").replace("]", ""),
                                                            str(sample.data.DP).replace(" ", "").replace("\'", "").replace("[", "").replace("]", ""),
                                                            str(sample.data.F1R2).replace(" ", "").replace("\'", "").replace("[", "").replace("]", ""), 
                                                            str(sample.data.F2R1).replace(" ", "").replace("\'", "").replace("[", "").replace("]", ""), 
                                                            str(sample.data.FAD).replace(" ", "").replace("\'", "").replace("[", "").replace("]", ""), 
                                                            str(sample.data.SB).replace(" ", "").replace("\'", "").replace("[", "").replace("]", "")  ] 
                                                            ) )
    return print_line


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
    os.system (command)
    print ("\t", command)

##############################################################################################################

if __name__ == "__main__":

    ti_df = pd.read_csv ( kwargs["TUMOR_INTERVAL"], sep = "\t", header = None)
    vcf_reader_tumor = vcf.Reader(open( kwargs ["TUMOR_MUTECT2_VCF"] + ".gz", "rb")  )
    vcf_reader_other = vcf.Reader(open( kwargs ["OTHER_MUTECT2_VCF"] + ".gz", "rb")  )

    print (kwargs)

    other_rescue_vcf = open ( kwargs ["OTHER_RESCUE_VCF"] , "a" )
    tumor_shared_variant = open ( kwargs ["TUMOR_SHARED_VARIANT_VCF"] , "a" )
    other_shared_variant = open ( kwargs ["OTHER_SHARED_VARIANT_VCF"] , "a" )
    tumor_unique_vcf = open ( kwargs ["TUMOR_UNIQUE_VCF"] , "a" )
    other_unique_vcf = open ( kwargs ["OTHER_UNIQUE_VCF"] , "a" )

    for k in range (ti_df.shape[0]):
        CHR, START, END, REF, ALT = ti_df.iloc[k][0], int(ti_df.iloc[k][1]), int(ti_df.iloc[k][2]), ti_df.iloc[k][3], ti_df.iloc[k][4]

        #1. bam file로 alt 가 1개 이상이라도 있는지 확인하기 (pysam)
        ALT_TYPE = "Ins" if len(ALT) > len(REF) else ("Del" if len(ALT) < len(REF) else ALT)
        if (len(ALT) == len(REF)) & (len(ALT) >= 2):  #이런건 지금 처리 못한다
            continue 
        Count = Variant ( kwargs["Sample_ID"], "Dura", kwargs["CASE_BAM_PATH"], CHR, END, 0, ALT_TYPE  )
        print ( "\n\n{}:{}\t{}>{}".format (Count.chr_, Count.pos_, REF, ALT), end = "\t")
        Count.pysam_read ( **kwargs )
        Count_blood = Variant ( kwargs["Sample_ID"], "Blood", kwargs["CONTROL_BAM_PATH"], CHR, END, 0, ALT_TYPE  )
        Count_blood.pysam_read ( **kwargs )


        #2. tumor mutect2 vcf에서 해당 부위 정보를 꺼내오기 (pyvcf)
        if (Count.base_counts [Count.alt_] >= kwargs["MINIMUM_ALT"]) :        # Dura에서는 일정 숫자 이상이고, Blood에서는 일정 숫자 이하일 것
            vcf_reader_tumor_fetch = vcf_reader_tumor.fetch ( CHR, START, END )
            
            for tumor_line in vcf_reader_tumor_fetch:
                # if "STR" in tumor_line.INFO.keys():   # Tumor VCF에서 STR==True라고 판명나면 그냥 버리기
                #     if tumor_line.INFO["STR"] == True:
                #         print (tumor_line.INFO["STR"])
                #         continue

                # BAMSNAP 찍기 (Dura와 Tumor 한번에)
                print ( "\t{} : {} / {}".format (kwargs ["TISSUE"], Count.base_counts [Count.alt_], Count.depth_ ) )
                print ( "\tBlood : {} / {}".format (Count_blood.base_counts [Count_blood.alt_], Count_blood.depth_ ) )
                if (Count_blood.base_counts [Count_blood.alt_] > kwargs["MINIMUM_ALT"]):
                    continue

                print ("\tBAMSNAP_IGV : {}".format ( kwargs["OTHER_SHARED_VARIANT_VCF"].replace ("vcf", str(CHR) + ":" + str(END) + ".jpg" ) ) )
                BAMSNAP_IGV ( BAM_DIR_LIST = " ".join ( [ kwargs["CASE_BAM_PATH"], kwargs["CASE_BAM_PATH"].replace( kwargs["TISSUE"], "Tumor"), kwargs["CONTROL_BAM_PATH"] ] ), TITLE_LIST = str(CHR) + ":" + str(END), POS = str(CHR) + ":" + str(END), OUT = kwargs["OTHER_SHARED_VARIANT_VCF"].replace ("vcf", str(CHR) + ":" + str(END) + ".jpg") )


                tumor_line.FILTER = "RESCUE"
                # print Tumor line
                print_tumor_line = print_line_vcf ( tumor_line )
                if ("Dura" in kwargs ["TISSUE"]) | ("Cortex" in kwargs ["TISSUE"]) | ("Ventricle" in kwargs ["TISSUE"]):
                    #print( "\t".join( [str(i) for i in print_tumor_line] ) )
                    print( "\t".join( [str(i) for i in print_tumor_line] ), file = tumor_shared_variant )

                # print other (Dura, Cortex, Ventricle) line
                other_line = tumor_line
                new_call_data = vcf.model.make_calldata_tuple(["GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"])(
                    GT = tumor_line.samples[1].data.GT,
                    AD = str( Count.depth_ -  Count.base_counts [Count.alt_] ) + "," + str( Count.base_counts [Count.alt_] ),
                    AF = round( Count.base_counts [Count.alt_] / Count.depth_ , 3),
                    DP = Count.depth_,  # Modify DP attribute
                    F1R2 = tumor_line.samples[1].data.F1R2,
                    F2R1 = tumor_line.samples[1].data.F2R1,
                    FAD = tumor_line.samples[1].data.FAD,
                    SB = tumor_line.samples[1].data.SB
                )
                other_line.samples[1].data = new_call_data
                print_rescue_line = print_line_vcf ( other_line )
                
                try:  
                    if  ( len ( list ( vcf_reader_other.fetch ( CHR, START, END ) ) ) == 0):      # 기존에 없어야 rescue로 새로 추가를 해주지
                        print( "\t".join( [str(i) for i in print_rescue_line] ), file = other_rescue_vcf )    
                        print( "\t".join( [str(i) for i in print_rescue_line] ), file = other_shared_variant )
                    else:       # 기존에 있다면 기존 것을 shared에만 출력해준다
                        for line_original in vcf_reader_other.fetch ( CHR, START, END ):
                            line_original.FILTER = "PASS"
                            print_line_original = print_line_vcf ( line_original )
                            print( "\t".join( [str(i) for i in print_line_original] ), file = other_shared_variant )

                except:  # 아예 fetch가 안되는 것은 해당 chromosome으로 된 contig가 없다는 얘기
                    print( "\t".join( [str(i) for i in print_rescue_line] ), file = other_rescue_vcf )
                    print( "\t".join( [str(i) for i in print_rescue_line] ), file = other_shared_variant )
                
                
                # print( "\t".join( [str(i) for i in print_line] ), file = tumor_unique_vcf )
                # print( "\t".join( [str(i) for i in print_line] ), file = other_unique_vcf )
        else:
            print ("\tNo shared mutation\n")

        

    other_rescue_vcf.close()
    tumor_shared_variant.close()
    other_shared_variant.close()
    tumor_unique_vcf.close()
    other_unique_vcf.close()


  