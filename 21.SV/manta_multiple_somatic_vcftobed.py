import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Here is usage direction.')
parser.add_argument('--INPUT_PATH', default="/data/project/Meningioma/21.SV/02.manta_multiple/02.PASS/02.Diploid/220930.Manta.Diploid.PASS.vcf")
parser.add_argument('--OUTPUT_PATH', default="/data/project/Meningioma/21.SV/02.manta_multiple/02.PASS/02.Diploid/220930.Manta.Diploid.PASS.chr.vcf")
parser.add_argument('--PANDAS_DIR', default="/data/project/Meningioma/21.SV/02.manta_multiple/03.pandas")
parser.add_argument('--ID', default="220930")
parser.add_argument('--CONTROL_SAMPLE', default="220930_Blood")
#warnings.simplefilter (action = 'ignore', category = FutureWarning)

args = parser.parse_args()
INPUT_PATH = args.INPUT_PATH
OUTPUT_PATH = args.OUTPUT_PATH
PANDAS_DIR = args.PANDAS_DIR
ID = args.ID
CONTROL_SAMPLE = args.CONTROL_SAMPLE

input_file = open (INPUT_PATH, "r")
output_file = open (OUTPUT_PATH, "w")

matrix = []
matrix2 = []
chrlist = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
df = pd.DataFrame (columns = ['chr', 'start', 'end', 'SVtype', 'sample', 'GT', 'SVlen', 'Qual', 'SVtool'])

def extract(line, info_dict):
    chr = line[0]
    sample = str(ID)
    SVtool = "Manta"
    Qual = int(line[-1].split(':')[-1].split(',')[0]) + int(line[-1].split(':')[-1].split(',')[1])

    if 'END' in info_dict:
        end = info_dict["END"]
    else:
        end = 0

    return [chr, 0, end, info_dict["SVTYPE"], sample, line[-1].split(':')[0] , 0, Qual, SVtool]


def inversion(line, info_dict):
    matrix2 = extract(line, info_dict)

    start = int(line[1])
    end = line[4].split(':')[1]

    if '[' in end:
        end = int(end.split('[')[0])
    elif ']' in end:
        end = int(end.split(']')[0])
    else:
        end = int(end)

    if end < start:
        start, end = end, start
    matrix2[1] = int(start); matrix2[2] = int(end); matrix2[6] = int(int(end) - int(start))

    return matrix2


def deldup(line, info_dict):
    matrix2 = extract(line, info_dict)

    start = int(line[1])
    end = int(matrix2[2])

    if end < start:
        start, end = end, start
    matrix2[1] = start; matrix2[2] = end; matrix2[6] = int(int(end) - int(start))

    return matrix2



sampledict_Ind_to_Name = {}
sampledict_Name_to_Ind = {}

while True:
    line = input_file.readline()
    line = line.rstrip('\n')

    if len(line) < 2:
        break

    if line.startswith("#") == True:                      # Header만 따로 저장이 필요할 경우
        print (line, file = output_file)
        if line.startswith("#CHROM") == True:
            line = line.split()
            for i in range (9, len(line)):
                sampledict_Ind_to_Name [i] = line[i]
                sampledict_Name_to_Ind[ line[i]  ] = i
            print (sampledict_Ind_to_Name)
        continue


    else:
        line = line.split()

        # 7
        info_list = line[7].split(';')
        if "IMPRECISE" in info_list:
            continue
        if (line[-1].split(':')[-1].split(',')[0] == 0) | (line[-1].split(':')[-1].split(',')[1] == 0):  ## SR이 하나롣 0이 있는경우 넘긴다
            continue
        info_dict = {}
        for i in range(len(info_list)):
            if '=' in info_list[i]:
                info_dict [info_list[i].split('=')[0]] = info_list[i].split('=')[1]

        #8. GT:FT:GQ:PL:PR:SR
        format_list = line[8].split(":")
        format_dict_Name_to_Ind, format_dict_Ind_to_Name = {}, {}
        for i in range(len( format_list )):
            format_dict_Name_to_Ind [ format_list[i] ] = i
            format_dict_Ind_to_Name [ i ] = format_list[i]

        CONTROL_GT =  line [ sampledict_Name_to_Ind [CONTROL_SAMPLE] ].split(":")[  format_dict_Name_to_Ind [ "GT"] ] 



        if CONTROL_GT in ["0/0"]:      # Blood에서 SV 없는 경우만 보고싶다
            if info_dict["SVTYPE"] == "BND":
                if line[0] in chrlist:
                    print (*line, sep = '\t', file = output_file)

                if line[0] in line[4]:         # 같은 chromosome인 경우
                    matrix2 = inversion(line, info_dict)
                else:
                    matrix2 = []               # translocation은 일단 넘기자
                    continue
            if (info_dict["SVTYPE"] == "DEL") | (info_dict["SVTYPE"] == "DUP"):
                matrix2 = deldup(line, info_dict)

            matrix.append(matrix2)




df = df.append(pd.DataFrame.from_records(matrix, columns = ['chr', 'start', 'end', 'SVtype', 'sample', 'GT', 'SVlen', 'Qual', 'SVtool']))


df = df[df["chr"].isin(chrlist)]
df["start"] = df["start"].astype(int)
df["end"] = df["end"].astype(int)
df["SVlen"] = df["SVlen"].astype(int)
df["Qual"] = df["Qual"].astype(int)
df= df.drop_duplicates()


print (df)

input_file.close()
output_file.close()


df_short_deletions = df[(df["SVtype"] == "DEL") & (df["SVlen"] < 10000)]
df_large_deletions = df[(df["SVtype"] == "DEL") & (df["SVlen"] >= 1000) & (df["SVlen"] <= 500000) ]
df_extra_deletions = df[(df["SVtype"] == "DEL") & (df["SVlen"] >= 300000)]
df_short_duplications = df[(df["SVtype"] == "DUP") & (df["SVlen"] < 10000)]
df_large_duplications = df[(df["SVtype"] == "DUP") & (df["SVlen"] >= 1000) & (df["SVlen"] <= 500000)]
df_extra_duplications = df[(df["SVtype"] == "DUP") & (df["SVlen"] >= 300000)]
df_short_inversions = df[(df["SVtype"] == "INV") & (df["SVlen"] < 10000)]
df_large_inversions = df[(df["SVtype"] == "INV") & (df["SVlen"] >= 1000) & (df["SVlen"] <= 500000)]

df_short_deletions.to_csv(PANDAS_DIR + "/"+ str(ID) + ".1.short_del.manta.bed", sep = '\t', index = False, header = False)
df_large_deletions.to_csv(PANDAS_DIR + "/"+ str(ID) + ".2.large_del.manta.bed", sep = '\t', index = False, header = False)
df_extra_deletions.to_csv(PANDAS_DIR + "/"+ str(ID) + ".3.extra_del.manta.bed", sep = '\t', index = False, header = False)
df_short_duplications.to_csv(PANDAS_DIR + "/"+ str(ID) + ".4.short_dup.manta.bed", sep = '\t', index = False, header = False)
df_large_duplications.to_csv(PANDAS_DIR + "/"+ str(ID) + ".5.large_dup.manta.bed", sep = '\t', index = False, header = False)
df_extra_duplications.to_csv(PANDAS_DIR + "/"+ str(ID) + ".6.extra_dup.manta.bed", sep = '\t', index = False, header = False)
