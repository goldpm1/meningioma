import argparse
import numpy as np
import pandas as pd
import random
import os
import gzip


parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--INPUT_VCF', type=str, default="/home/goldpm1/Meningioma/04.mutect/02.PASS/221026_Dura.MT2.FMC.vcf")
parser.add_argument('--OUTPUT_VCF', type=str, default="/home/goldpm1/Meningioma/04.mutect/02.PASS/221026_Dura.MT2.FMC.HF.vcf")
parser.add_argument('--SAMPLE_THRESHOLD', type=str, default="Dura,Tumor")
parser.add_argument('--DP_THRESHOLD', type=int, default=30)
parser.add_argument('--ALT_THRESHOLD', type=int, default=1)
parser.add_argument('--REMOVE_MULTIALLELIC', type=bool, default=True)
parser.add_argument('--PASS', type=bool, default=True)
parser.add_argument('--REMOVE_MITOCHONDRIAL_DNA', type=bool, default=True)

args = parser.parse_args()

kwargs = {}
kwargs["INPUT_VCF"] = args.INPUT_VCF
kwargs["OUTPUT_VCF"] = args.OUTPUT_VCF
kwargs["SAMPLE_THRESHOLD"] = args.SAMPLE_THRESHOLD.split(",")
kwargs["DP_THRESHOLD"] = args.DP_THRESHOLD
kwargs["ALT_THRESHOLD"] = args.ALT_THRESHOLD
kwargs["REMOVE_MULTIALLELIC"] = args.REMOVE_MULTIALLELIC
kwargs["PASS"] = args.PASS
kwargs["REMOVE_MITOCHONDRIAL_DNA"] = args.REMOVE_MITOCHONDRIAL_DNA



def parsing (line):
    CHR, POS, REF, ALT = line[0], int(line[1]), line[3], line[4]

    # 7. info
    info_list = line[7].split(';')
    info_dict = {}
    for i in range( len(info_list) ):
        if "=" in info_list[i]:
            info_dict [info_list[i].split('=')[0]] = info_list[i].split('=')[1]

    # 8 : format
    format_list = line[8].split(":")

    # 9 ~ : sample
    sample_dict = []
    for i in range(9, len(line)):
        temp_dict = {}
        if line[i] == ".":      # 없는 경우
            for j_index, j in enumerate (format_list):
                temp_dict[ format_list[j_index] ] = "."
            sample_dict.append (temp_dict)
        else:
            for j_index, j in enumerate(line[i].split(":")):
                temp_dict[ format_list[j_index] ] = j
            sample_dict.append (temp_dict)

    return CHR, POS, REF, ALT, info_list, info_dict, format_list, sample_dict



input_file = open ( kwargs["INPUT_VCF"], "r")
output_file = open ( kwargs["OUTPUT_VCF"], "w")
sample_name = []

for line in input_file.readlines():
    if "gz" in kwargs["INPUT_VCF"]:
        line = line.decode('utf-8')
    line = line.rstrip("\n")
    
    if line[0:4] == "#CHR":
        sample_name = line.split("\t")[9:]     # multisample일때 samplename을 list 화
    if line[0] == "#": # Header 저장
        print (line, file = output_file) # header 그대로 출력
        continue

    line = line.split("\t")

    CHR, POS, REF, ALT , info_list, info_dict, format_list, sample_dict = parsing(line)
    ## info_list : ['AS_FilterStatus=base_qual|weak_evidence,base_qual', 'AS_SB_TABLE=27,17|10,11|4,3', 'DP=88', 'ECNT=1', 'GERMQ=1', 'MBQ=20,20,20', 'MFRL=177,175,205', 'MMQ=60,60,60', 'MPOS=50,28', 'NALOD=-1.054e+00,0.623', 'NLOD=0.188,3.94', 'PON', 'POPAF=0.613,1.39', 'RPA=17,16,18', 'RU=A', 'STR', 'STRQ=1', 'TLOD=11.20,5.37']
    ## info_dict : {'AS_FilterStatus': 'base_qual|weak_evidence,base_qual', 'AS_SB_TABLE': '27,17|10,11|4,3', 'DP': '88', 'ECNT': '1', 'GERMQ': '1', 'MBQ': '20,20,20', 'MFRL': '177,175,205', 'MMQ': '60,60,60', 'MPOS': '50,28', 'NALOD': '-1.054e+00,0.623', 'NLOD': '0.188,3.94', 'POPAF': '0.613,1.39', 'RPA': '17,16,18', 'RU': 'A', 'STRQ': '1', 'TLOD': '11.20,5.37'}
    ## format_list : ['GT', 'AD', 'AF', 'DP', 'F1R2', 'F2R1', 'FAD', 'SB']
    ## sample_dict [] :  [{'GT': '0|0', 'AD': '63,0', 'AF': '0.022', 'DP': '63', 'F1R2': '16,0', 'F2R1': '23,0', 'FAD': '44,0', 'PGT': '0|1', 'PID': '89561722_C_G', 'PS': '89561722', 'SB': '30,33,0,0'}, {'GT': '0|1', 'AD': '126,5', 'AF': '0.041', 'DP': '131', 'F1R2': '37,0', 'F2R1': '47,2', 'FAD': '93,3', 'PGT': '0|1', 'PID': '89561722_C_G', 'PS': '89561722', 'SB': '56,70,2,3'}]
    ## sample_name[] :  ['221026_Blood', '221026_Dura']

    

    check = 0
    if kwargs["REMOVE_MITOCHONDRIAL_DNA"] == True:
        if CHR not in ["chr" + str(i) for i in range (1,23) ] + ["chrX, chrY"]:
            check = 1        
    if kwargs["PASS"] == True:
        if "PASS" not in line[6]:
            check = 1
    if kwargs["REMOVE_MULTIALLELIC"]  == True:
        if "," in line[4]:
            check = 1
        
    for i in range (9, len (line)):
        #print ("sample = {}  → DP = {}, AD = {}".format (sample_name[i - 9], sample_dict [i - 9] ["DP"], sample_dict [i - 9] ["AD"]))
        check2 = False
        if kwargs["SAMPLE_THRESHOLD"] == "all":
            check2 = True
        else:
            for SAMPLE_THESHOLD_samplename in kwargs["SAMPLE_THRESHOLD"]:  # [Dura, Tumor] 에 포함되는 sample만 DP, AD threshold를 판단한다
                if SAMPLE_THESHOLD_samplename in sample_name[i - 9]:
                    check2 = True
        
        if check2 == True:
            if int( sample_dict [i - 9] ["DP"] ) < kwargs["DP_THRESHOLD"]:
                check = 1
                break
            else:
                if ( int ( sample_dict [i - 9] ["AD"].split(",")[0] ) < kwargs ["ALT_THRESHOLD"] ) | ( int ( sample_dict [i - 9] ["AD"].split(",")[1] ) < kwargs ["ALT_THRESHOLD"] ):
                    check = 1
                    break
    
    if check == 0:       # 조건에 안 걸릴 때에만 출력
        print ("\t".join( line ), file = output_file)       
        
input_file.close()
output_file.close()
            
    
