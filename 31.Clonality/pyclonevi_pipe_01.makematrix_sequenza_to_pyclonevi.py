import pandas as pd
import pybedtools
import pysam
import argparse
import random

parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--Sample_ID', type=str, default="220930")
parser.add_argument('--TISSUE', type=str, default="Tumor")
parser.add_argument('--SEQUENZA_MUTATION_PATH', type=str, default="/home/goldpm1/Meningioma/11.cnv/2.sequenza/220930_Tumor_mutations.txt")
parser.add_argument('--SEQUENZA_SEGMENT_PATH', type=str, default="/home/goldpm1/Meningioma/11.cnv/2.sequenza/220930_Tumor_segments.txt")
parser.add_argument('--SEQUENZA_TO_PYCLONEVI_MATRIX_PATH', type=str, default="/home/goldpm1/Meningioma/31.Clonality/01.make_matrix/220930/220930.sequenza_to_pyclonevi.tsv")
parser.add_argument('--SEQUENZA_PLOIDY_PATH', type=str, default="/home/goldpm1/Meningioma/11.cnv/2.sequenza/220930_Tumor_confints_CP.txt")
parser.add_argument('--SEQUENZA_PURITY_PLOIDY_PATH', type=str, default="/home/goldpm1/Meningioma/11.cnv/2.sequenza/220930_Tumor_confints_CP.txt")
parser.add_argument('--MUTECT_OUTPUT_PATH', type=str, default="/home/goldpm1/Meningioma/04.mutect/02.PASS/220930_Tumor.MT2.FMC.HF.vcf")
parser.add_argument('--HC_OUTPUT_PATH', type=str, default="/home/goldpm1/Meningioma/06.hc/03.HF/220930/Blood/220930_Blood.DP100.vcf")
parser.add_argument('--HC_BLOOD_RANDOM_PICK_PATH', type=str, default="/home/goldpm1/Meningioma/31.Clonality/01.make_matrix/220930/220930.HC.random_pick_50.bed")
parser.add_argument('--GVCF_OUTPUT_PATH', type=str, default="/home/goldpm1/Meningioma/05.gvcf/02.remove_nonref/220930/220930_Tumor.g.vcf.gz")

args = parser.parse_args()

SEQUENZA_MUTATION_PATH = args.SEQUENZA_MUTATION_PATH
SEQUENZA_SEGMENT_PATH = args.SEQUENZA_SEGMENT_PATH
SEQUENZA_PLOIDY_PATH = args.SEQUENZA_PLOIDY_PATH
SEQUENZA_PURITY_PLOIDY_PATH=args.SEQUENZA_PURITY_PLOIDY_PATH
SEQUENZA_TO_PYCLONEVI_MATRIX_PATH = args.SEQUENZA_TO_PYCLONEVI_MATRIX_PATH
MUTECT_OUTPUT_PATH = args.MUTECT_OUTPUT_PATH
GVCF_OUTPUT_PATH = args.GVCF_OUTPUT_PATH
HC_OUTPUT_PATH = args.HC_OUTPUT_PATH
HC_BLOOD_RANDOM_PICK_PATH = args.HC_BLOOD_RANDOM_PICK_PATH
TISSUE = args.TISSUE
Sample_ID = args.Sample_ID

tbx = pysam.TabixFile (GVCF_OUTPUT_PATH)

SAMPLE_ID = Sample_ID + "_" + TISSUE

ploidy_df = pd.read_csv (SEQUENZA_PLOIDY_PATH, sep = "\t")
TUMOR_PURITY = ploidy_df.iloc[0]["cellularity"]
TOTAL_PLOIDY=ploidy_df.iloc[0]["ploidy.estimate"]

mutation_df = pd.read_csv (SEQUENZA_MUTATION_PATH, sep = "\t")  
segment_df = pd.read_csv (SEQUENZA_SEGMENT_PATH, sep = "\t")

mutation_df["end"] = mutation_df["position"]
mutation_df["start"] = mutation_df["position"] - 1
if "map_ratio" in mutation_df.columns():  #hgt19to38의 경우 map_ratio라는 쓸데없는 column이 있다
    mutation_df = mutation_df.drop ('map_ratio', axis = 1)

 
mutation_df = mutation_df.reindex (columns = ["chromosome", "start", "end", "GC.percent", "good.reads", "adjusted.ratio", "F", "mutation", "CNn", "CNt", "Mt"] )
mutation_df.head()


######## Mutect call을 골라줌 ###################
mutect_pybed_object = pybedtools.BedTool(MUTECT_OUTPUT_PATH)
segment_pybed_object =  pybedtools.BedTool.from_dataframe(segment_df)

ab = segment_pybed_object.intersect(mutect_pybed_object, wb = True)
for index, contents in enumerate( ab[0] ):
    if "GT:" in contents:
        parsing_sample_index = index 
        break


output_file = open(SEQUENZA_TO_PYCLONEVI_MATRIX_PATH, "a")

for interval in ab:
    CHR, POS = interval[0], interval[2]
    line_dict = {}

    line_dict ["mutation_id"] = str(CHR) + "_"  + str(POS)
    line_dict ["sample_id"] = SAMPLE_ID
    line_dict ["ref_counts"] = str( interval[parsing_sample_index + 2].split(":")[1].split(",")[0])
    line_dict ["alt_counts"] = str( interval[parsing_sample_index + 2].split(":")[1].split(",")[1])
    line_dict["normal_cn"] = str ( 2 )
    line_dict["major_cn"] = str ( interval[10] )
    line_dict["minor_cn"] = str ( interval[11] )
    line_dict["tumour_content"] = str( TUMOR_PURITY  )
    line_dict["type"] = "Mutect2"

    #print ( interval [parsing_sample_index - 1 : parsing_sample_index + 2])
    line_dict["gene"] = str( interval [parsing_sample_index - 1] ).split(";")[-1].split("|")[3]
    line_dict["Variant_Classification"] = str( interval [parsing_sample_index - 1] ).split(";")[-1].split("|")[1]
    
    #print ( '\t'.join ( list( line_dict.values()) ) )
    print ( '\t'.join ( list( line_dict.values()) ), file = output_file )
        
    
    
# ######## HC call을 골라줌 ###################
# hc_blood_pybed_object = pybedtools.BedTool ( HC_BLOOD_RANDOM_PICK_PATH )   # BED
# hc_pybed_object = pybedtools.BedTool(HC_OUTPUT_PATH)  # VCF

# hc_selected_pybed_object = hc_blood_pybed_object.intersect(hc_pybed_object, wa = True, wb = True)       # Dura/Tumor에서 발견된 HC call 중에서 선택된 50개만 고름 (concordance가 95%이기 때문에 47개정도 예상)

# ab = hc_selected_pybed_object.intersect(segment_pybed_object, wa = True, wb = True)

# for index, contents in enumerate( ab[0] ):      # GT:AD 다음 칸의 정보를 얻는다
#     if "GT:" in contents:
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
#     line_dict["major_cn"] = str ( interval[ -3 ] )
#     line_dict["minor_cn"] = str ( interval[ -2 ] )
#     line_dict["tumour_content"] = str( TUMOR_PURITY  )
#     line_dict["type"] = "HC"
#     line_dict["gene"] = str( interval [parsing_sample_index - 1] ).split(";")[-1].split("|")[3]
#     line_dict["Variant_Classification"] = str( interval [parsing_sample_index - 1] ).split(";")[-1].split("|")[1]
    
#     #print ( '\t'.join ( list( line_dict.values()) ) )
#     print ( '\t'.join ( list( line_dict.values()) ), file = output_file )
   

output_file.close()

pd.DataFrame ( { "Purity" : [ TUMOR_PURITY ], " Ploidy" : [ TOTAL_PLOIDY ] }).to_csv ( SEQUENZA_PURITY_PLOIDY_PATH, sep = "\t", index = False, header = True)