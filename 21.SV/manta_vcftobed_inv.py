import pandas as pd
import numpy as np
import csv
import gzip
import glob
import time
import sys

InputPath = sys.argv[1]
PANDAS_DIR = sys.argv[2]
ID = sys.argv[3]


input_file = open (InputPath, "r")


matrix = []
matrix2 = []
chrlist = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
df = pd.DataFrame (columns = ['chr', 'start', 'end', 'SVtype', 'sample', 'GT', 'SVlen', 'Qual', 'SVtool'])

def extract(line, info_dict):
    chr = line[0]
    sample = str (ID)
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
    end = int(matrix2[2])

    if end < start:
        start, end = end, start
    matrix2[1] = int(start); matrix2[2] = int(end); matrix2[6] = int(int(end) - int(start))

    return matrix2





while True:
    line = input_file.readline()
    line = line.rstrip('\n')

    if len(line) < 2:
        break

    if line[0] == "#":                      # Header만 따로 저장이 필요할 경우
        continue


    line = line.split()

    info_list = line[7].split(';')
    if "IMPRECISE" in info_list:
        continue
    if (line[-1].split(':')[-1].split(',')[0] == 0) | (line[-1].split(':')[-1].split(',')[1] == 0):  ## SR이 하나롣 0이 있는경우 넘긴다
        continue

    info_dict = {}
    for i in range(len(info_list)):
        if '=' in info_list[i]:
            info_dict [info_list[i].split('=')[0]] = info_list[i].split('=')[1]


    if info_dict["SVTYPE"] == "INV":      # 일단 translocation은 넘기자
        matrix2 = inversion(line, info_dict)


    matrix.append(matrix2)




df = df.append(pd.DataFrame.from_records(matrix, columns = ['chr', 'start', 'end', 'SVtype', 'sample', 'GT', 'SVlen', 'Qual', 'SVtool']))


df = df[df["chr"].isin(chrlist)]
df["start"] = df["start"].astype(int)
df["end"] = df["end"].astype(int)
df["SVlen"] = df["SVlen"].astype(int)
df["Qual"] = df["Qual"].astype(int)
df= df.drop_duplicates()

df_inv = df[(df["SVtype"] == "INV")]


#print (df_short_deletions)

input_file.close()

df_short_inversions = df_inv[(df_inv["SVlen"] < 10000)]
df_large_inversions = df_inv[(df_inv["SVlen"] >= 1000)]

df_short_inversions.to_csv(PANDAS_DIR + "/" + str(ID) + ".7.short_inv.manta.bed", sep = '\t', index = False, header = False)
df_large_inversions.to_csv(PANDAS_DIR + "/" + str(ID) + ".8.large_inv.manta.bed", sep = '\t', index = False, header = False)
