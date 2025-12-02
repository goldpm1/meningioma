import pandas as pd
import gzip, os, vcf
import argparse
import pybedtools
import warnings
warnings.simplefilter (action = 'ignore')
import re

# 숫자 추출 함수 정의
def extract_chr_num(chr_str):
    match = re.search(r'chr(\d+|X|Y|M)', chr_str)
    if match:
        val = match.group(1)
        # X, Y, M 처리
        if val == 'X':
            return 23
        elif val == 'Y':
            return 24
        elif val in ('M', 'MT'):
            return 25
        else:
            return int(val)
    else:
        return float('inf')  # 알 수 없는 경우 맨 뒤로
    



parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--Sample_ID', type=str, default="230405_2")
parser.add_argument('--TISSUE', type=str, default="Tumor")
parser.add_argument('--FACETCNV_OUTPUT_PATH', type=str, default="/home/goldpm1/Meningioma/11.cnv/5.facetcnv/230405_2/Tumor/230405_2.vcf.gz")
parser.add_argument('--FACETCNV_TO_BED_DF_PATH', type=str, default="")
parser.add_argument('--FACETCNV_TO_PYCLONEVI_MATRIX_PATH', type=str, default="/home/goldpm1/Meningioma/31.Clonality/01.make_matrix/230405_2/230405_2.facetcnv_to_pyclonevi.tsv")
parser.add_argument('--FACETCNV_PURITY_PLODY_PATH', type=str, default="/home/goldpm1/Meningioma/31.Clonality/01.make_matrix/230405_2/230405_2.facetcnv_to_pyclonevi.tsv")
parser.add_argument('--MUTECT_PATH', type=str, default="/home/goldpm1/Meningioma/04.mutect/02.PASS/230405_2_Tumor.MT2.FMC.HF.vcf")
parser.add_argument('--REMOVE_TEMP', type=str, default="True")
# parser.add_argument('--HC_OUTPUT_PATH', type=str, default="/home/goldpm1/Meningioma/06.hc/01.call/230405_2/Tumor/230405_2_Tumor.vcf")
# parser.add_argument('--HC_BLOOD_RANDOM_PICK_PATH', type=str, default="/home/goldpm1/Meningioma/31.Clonality/01.make_matrix/230405_2/230405_2.HC.random_pick_50.bed")

args = parser.parse_args()

kwargs = {}
Sample_ID = args.Sample_ID
TISSUE = args.TISSUE
FACETCNV_OUTPUT_PATH = args.FACETCNV_OUTPUT_PATH
FACETCNV_TO_BED_DF_PATH=args.FACETCNV_TO_BED_DF_PATH
FACETCNV_PURITY_PLODY_PATH=args.FACETCNV_PURITY_PLODY_PATH
MUTECT_PATH = args.MUTECT_PATH
# HC_OUTPUT_PATH=args.HC_OUTPUT_PATH
# HC_BLOOD_RANDOM_PICK_PATH=args.HC_BLOOD_RANDOM_PICK_PATH
FACETCNV_TO_PYCLONEVI_MATRIX_PATH = args.FACETCNV_TO_PYCLONEVI_MATRIX_PATH
REMOVE_TEMP = args.REMOVE_TEMP

# GVCF_OUTPUT_PATH="/home/goldpm1/Meningioma/05.gvcf/01.call/221026/221026_Tumor.g.vcf.gz"
# tbx = pysam.TabixFile (GVCF_OUTPUT_PATH)

SAMPLE_ID = Sample_ID + "_" + TISSUE



################## FACETCNV가 준 vcf.gz 가지고 Segment bed object 만들기 ##############

if "gz" in FACETCNV_OUTPUT_PATH:
    input_file = gzip.open (FACETCNV_OUTPUT_PATH, "r")
else:
    input_file = open (FACETCNV_OUTPUT_PATH, "r")
sample_name = []


def parsing (line):
    CHR, POS, REF, ALT = line[0], int(line[1]), line[3], line[4]
    # 7. info
    info_list = line[7].split(';')
    info_dict = {}
    for i in range( len(info_list) ):
        if "=" in info_list[i]:
            info_dict [ info_list[i].split('=')[0] ] = info_list[i].split('=')[1]

    return CHR, POS, REF, ALT, info_list, info_dict



# FACETCNV 의 vcf.gz를 읽어서 bed file로 바꾸는 것. 그런데 빈 공간이 있어서 주의

matrix = [] 
colnames = ['CHR','START','END', 'MAJOR_CN', 'MINOR_CN', 'NORMAL_CN', 'TUMOR_PURITY' ]
bed_df = pd.DataFrame (columns = colnames) 

for line in input_file.readlines():
    if "gz" in FACETCNV_OUTPUT_PATH:
	    line = line.decode('utf-8')
    line = line.rstrip('\n')

     
    
    if line[0] == "#": # Header 저장
        if line[0:8] == "##purity":
            TUMOR_PURITY = line.split("=")[1]
            if TUMOR_PURITY == "NA":
                TUMOR_PURITY = 1.0
            else:
                TUMOR_PURITY = round ( float(TUMOR_PURITY), 2 )
        elif line[0:8] == "##ploidy":
            TOTAL_PLOIDY = line.split("=")[1]
        continue
    
    else:
        line = line.split("\t")
        CHR, POS, REF, ALT , info_list, info_dict = parsing(line)
        START, END = POS, info_dict["END"]
        try:   #"LCN_EM에 . 같은 게 있어서 에러날 때가 있다"
            MAJOR_CN, MINOR_CN =  int(info_dict["TCN_EM"])  - int (info_dict["LCN_EM"]),  int (info_dict["LCN_EM"])
            NORMAL_CN = 2           
            if (MAJOR_CN == 0) & (MINOR_CN == 0): #230419_Tumor는 chr22에서 이상하게 TCN_EM=0;LCN_EM=0;로 나와서 pyclone-vi를 못 돌리게 한다
                MAJOR_CN, MINOR_CNt = 1, 0
        except:
            continue
        output_line = [str(CHR), str(START), str(END), str(MAJOR_CN), str(MINOR_CN), str(NORMAL_CN),  str(TUMOR_PURITY)  ]
        output_dict = {}
        for colname_index, colname in enumerate (colnames):
            output_dict[colname] = output_line[colname_index]
            
        bed_df = bed_df.append ( pd.Series( output_dict ), ignore_index = True )
        #bed_df = pd.concat ( [bed_df,  pd.Series( output_dict )], axis = 0)


################### Gap interval도 채워서 넣어주기 #####################################
# 이전 end 값 계산
bed_df ["prev_END"] = bed_df .groupby("CHR")["END"].shift(1)    # 즉, 같은 chromosome 안에서 바로 이전 구간의 END 위치를 현재 행에 기록함.
bed_df ["prev_CHR"] = bed_df ["CHR"].shift(1)

# gap이 있는 경우만 필터링
gaps = bed_df [(bed_df ["CHR"] == bed_df ["prev_CHR"]) & (bed_df ["START"] > bed_df ["prev_END"]) ].copy()

# gap interval 만들기
gap_bed = pd.DataFrame({
    "CHR": gaps["CHR"],    "START": gaps["prev_END"],     "END": gaps["START"],
    "MAJOR_CN": 1,    "MINOR_CN": 1,    "NORMAL_CN": 2,     "TUMOR_PURITY": 1.0
})

bed_df = bed_df.drop ( columns = ["prev_END", "prev_CHR" ] ) # 원래 bed_f에서 보조 컬럼 제거
merged_df = pd.concat ( [ bed_df, gap_bed ], ignore_index = True) # 기존 bed_f와 gap df 합치기


# 자연스러운 chr 순서로 정렬
merged_df = merged_df.astype ( {'START':'int' , 'END' : 'int' , 'MAJOR_CN':'int', "MINOR_CN" : "int", "NORMAL_CN" : "int", "TUMOR_PURITY" : "float" } )
merged_df["CHR_sort"] = merged_df["CHR"].apply(extract_chr_num)
merged_df = merged_df.sort_values(by=["CHR_sort", "START"]).drop(columns="CHR_sort").reset_index(drop=True)

bed_df = merged_df





bed_df.to_csv (FACETCNV_TO_PYCLONEVI_MATRIX_PATH + ".temp1", sep = "\t", index = False, header = False)
bed_df.to_csv (FACETCNV_TO_BED_DF_PATH, sep = "\t", index = False, header = True)
pd.DataFrame ( { "Purity" : [ TUMOR_PURITY ], " Ploidy" : [ TOTAL_PLOIDY ] }).to_csv ( FACETCNV_PURITY_PLODY_PATH, sep = "\t", index = False, header = True)






####################################################################################################
######## Mutect call을 골라줌 ###################
bed_pybed_object = pybedtools.BedTool.from_dataframe(bed_df)                   # CNV 정보를 가지고 있는 segment 정보
# mutect_pybed_object = pybedtools.BedTool( MUTECT_PATH )       # mutect으로 call한 mutation
# a = bed_pybed_object.intersect(mutect_pybed_object, wb = True) 


os.system ("bedtools intersect -a " + MUTECT_PATH + " -b " + FACETCNV_TO_PYCLONEVI_MATRIX_PATH + ".temp1" +  " -wa -wb > " + FACETCNV_TO_PYCLONEVI_MATRIX_PATH + ".temp2")

vcf_reader = vcf.Reader(open (MUTECT_PATH , "r"))
for sampleindex, samplename in enumerate(vcf_reader.samples):
    if "Blood" not in samplename:
        #print ( samplename, sampleindex)       # Tumor, Dura, Falx가  vcf 파일에서 0번째인지 1번째인지 보는 것
        break



input_file = open (FACETCNV_TO_PYCLONEVI_MATRIX_PATH + ".temp2", "r")
output_file = open(FACETCNV_TO_PYCLONEVI_MATRIX_PATH, "a")

for interval in input_file.readlines():
    interval = interval.rstrip("\n").split("\t")

    
    line_dict = {}
    
    for index, contents in enumerate( interval ):
        if "GT:" in contents:
            parsing_sample_index = index 

            break

    line_dict ["mutation_id"] = str(interval [0]) + "_" + str(interval[1]) + "_" + str(interval[3]) + "_" + str(interval[4])
    line_dict ["sample_id"] = SAMPLE_ID
    line_dict ["ref_counts"] = str( interval[ parsing_sample_index + 1 + sampleindex ].split(":")[1].split(",")[0])     # GT 다음칸인지, 다다음칸인지는 sample name이 몇 번째에 오는지를 보고 판단해야 한다
    line_dict ["alt_counts"] = str( interval[ parsing_sample_index + 1 + sampleindex ].split(":")[1].split(",")[1])
    line_dict["normal_cn"] = str (interval[-2])
    line_dict["major_cn"] = str (interval[-4])
    line_dict["minor_cn"] = str (interval[-3])
    line_dict["tumour_content"] = str( interval[-1])
    line_dict["type"] = "Mutect2"
    line_dict["gene"] = str( interval [parsing_sample_index - 1] ).split(";")[-1].split("|")[3]
    line_dict["Variant_Classification"] = str( interval [parsing_sample_index - 1] ).split(";")[-1].split("|")[1]
    
    #print ( '\t'.join ( list( line_dict.values()) ) )
    print ( '\t'.join ( list( line_dict.values()) ), file = output_file )

input_file.close()

if REMOVE_TEMP == "True":
    os.system ("rm -rf " + FACETCNV_TO_PYCLONEVI_MATRIX_PATH + ".txt")
    os.system ("rm -rf " + FACETCNV_TO_PYCLONEVI_MATRIX_PATH + ".temp1")
    os.system ("rm -rf " + FACETCNV_TO_PYCLONEVI_MATRIX_PATH + ".temp2")
    



# ####################################################################################################
#  ######## HC call을 골라줌 ###################

# hc_blood_pybed_object = pybedtools.BedTool ( HC_BLOOD_RANDOM_PICK_PATH )  
# hc_pybed_object = pybedtools.BedTool(HC_OUTPUT_PATH)  #VCF

# hc_selected_pybed_object = hc_blood_pybed_object.intersect(hc_pybed_object, wa = True, wb = True)       # Dura/Tumor에서 발견된 HC call 중에서 선택된 50개만 고름 (concordance가 95%이기 때문에 47개정도 예상)


# ab = hc_selected_pybed_object.intersect(bed_pybed_object, wa = True, wb = True)


# for index, contents in enumerate( ab[0] ):
#     if "GT" in contents:
#         parsing_sample_index = index 
#         break


# for interval in ab:
#     CHR, POS = interval[0], int (interval[2]) - 1
#     line_dict = {}

#     line_dict ["mutation_id"] = str(CHR) + "_"  + str(POS)
#     line_dict ["sample_id"] = SAMPLE_ID
#     line_dict ["ref_counts"] = str( interval [parsing_sample_index + 1] ).split(":")[1].split(",")[0]
#     line_dict ["alt_counts"] = str( interval [parsing_sample_index + 1] ).split(":")[1].split(",")[1]
#     line_dict["normal_cn"] = str ( 2 )
#     line_dict["major_cn"] = str ( interval[ -4 ] )
#     line_dict["minor_cn"] = str ( interval[ -3 ] )
#     line_dict["tumour_content"] = str( TUMOR_PURITY  )
#     line_dict["type"] = "HC"
#     line_dict["gene"] = str( interval [parsing_sample_index - 1] ).split(";")[-1].split("|")[3]
#     line_dict["Variant_Classification"] = str( interval [parsing_sample_index - 1] ).split(";")[-1].split("|")[1]
    
#     #print ( '\t'.join ( list( line_dict.values()) ) )
#     print ( '\t'.join ( list( line_dict.values()) ), file = output_file )
   

output_file.close()
