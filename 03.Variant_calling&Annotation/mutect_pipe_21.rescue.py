import argparse
import numpy as np
import pandas as pd
import os, pysam, vcf, pybedtools
import gzip


parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--Sample_ID', type=str)
parser.add_argument('--TISSUE', type=str)
parser.add_argument('--NUM', type=str)
parser.add_argument('--MULTIPLE_VCF_PATH', type=str)
parser.add_argument('--MULTIPLE_VCF_GZ_PATH', type=str)
parser.add_argument('--INDIVIDUAL_VCF_PATH', type=str)
parser.add_argument('--INDIVIDUAL_VCF_GZ_PATH', type=str)
parser.add_argument('--INDIVIDUAL_RESCUED_VCF_PATH', type=str)
parser.add_argument('--INDIVIDUAL_RESCUED_VCF_GZ_PATH', type=str)
parser.add_argument('--INDIVIDUAL_UNIQUE_VCF_PATH', type=str)
parser.add_argument('--INDIVIDUAL_UNIQUE_VCF_GZ_PATH', type=str)

args = parser.parse_args()

kwargs = {}
kwargs["Sample_ID"] = args.Sample_ID
kwargs["TISSUE"] = args.TISSUE
kwargs["NUM"] = args.NUM
kwargs["MULTIPLE_VCF_PATH"] = args.MULTIPLE_VCF_PATH
kwargs["MULTIPLE_VCF_GZ_PATH"] = args.MULTIPLE_VCF_GZ_PATH
kwargs["INDIVIDUAL_VCF_PATH"] = args.INDIVIDUAL_VCF_PATH
kwargs["INDIVIDUAL_VCF_GZ_PATH"] = args.INDIVIDUAL_VCF_GZ_PATH
kwargs["INDIVIDUAL_RESCUED_VCF_PATH"] = args.INDIVIDUAL_RESCUED_VCF_PATH
kwargs["INDIVIDUAL_RESCUED_VCF_GZ_PATH"] = args.INDIVIDUAL_RESCUED_VCF_GZ_PATH
kwargs["INDIVIDUAL_UNIQUE_VCF_PATH"] = args.INDIVIDUAL_UNIQUE_VCF_PATH
kwargs["INDIVIDUAL_UNIQUE_VCF_GZ_PATH"] = args.INDIVIDUAL_UNIQUE_VCF_GZ_PATH




if __name__ == "__main__":
    pybed_individual = pybedtools.BedTool( kwargs["INDIVIDUAL_VCF_PATH"] )
    pybed_multiple = pybedtools.BedTool( kwargs["MULTIPLE_VCF_PATH"] )
    
    input_file_individual_vcf = vcf.Reader (open ( kwargs["INDIVIDUAL_VCF_PATH"], "r"))
    input_file_multiple_vcf = vcf.Reader (open ( kwargs["MULTIPLE_VCF_PATH"], "r"))
    sample_name_individual_vcf = input_file_individual_vcf.samples  # ['220930_Blood', '220930_Dura']
    sample_name_multiple_vcf = input_file_multiple_vcf.samples       # ['220930_Blood', '220930_Dura', '220930_Tumor']
    
    
    # Header는 individual 그대로 출력
    output_file_exists =  os.path.exists ( kwargs["INDIVIDUAL_RESCUED_VCF_PATH"] + ".unsorted"  )
    output_file_unique_exists =  os.path.exists ( kwargs["INDIVIDUAL_UNIQUE_VCF_PATH"] + ".unsorted" )

    print ("output_file_exists = {}".format(output_file_exists))
    print ("output_file_unique_exists = {}\n".format(output_file_unique_exists))

    if output_file_exists == False:
        output_file = open ( kwargs["INDIVIDUAL_RESCUED_VCF_PATH"] + ".unsorted" , "w" )
    else:
        output_file = open ( kwargs["INDIVIDUAL_RESCUED_VCF_PATH"] + ".unsorted" , "a" )   # 덧붙여주기

    if output_file_unique_exists == False:
        output_file_unique = open ( kwargs["INDIVIDUAL_UNIQUE_VCF_PATH"] + ".unsorted" , "w" )
    else:
        output_file_unique = open ( kwargs["INDIVIDUAL_UNIQUE_VCF_PATH"]  + ".unsorted" , "a" )    # 덧붙여주기
    
    input_file = open ( kwargs["INDIVIDUAL_VCF_PATH"], "r")
    for line in input_file.readlines():
        line = line.rstrip("\n")            
        if line [0] == '#':        # Header 저장
            if output_file_exists == False:        # 처음에 파일 없고 새로 생성됐을 때에만 Header 저장
                print (line, file = output_file)
            if output_file_unique_exists == False:
                print (line, file = output_file_unique)
            continue
        else:
            break
    input_file.close()
    
    
    
    # Mutect call ∩ Individual call  → 그냥 Individual call을 출력
    ab = pybed_individual.intersect(pybed_multiple)
    for line in ab:
        print ( "Mutect call ∩ Individual call", line[0:5], sep = "\t")
        print (*line, sep = "\t", file = output_file)
        
    # Individual call - Multiple call  → 그냥 Individual call을 출력. unique는 이것만 받는다
    ab = pybed_individual.intersect(pybed_multiple, v = True)
    for line in ab:
        print ( "Individual call - Multiple call", line[0:5], sep = "\t")
        print (*line, sep = "\t", file = output_file)
        print (*line, sep = "\t", file = output_file_unique)
        
        
    # Multiple call - Individual call → Multiple call을 이식해서 rescuedㅔ 출력해줌
    ab = pybed_multiple.intersect(pybed_individual, v = True)
    s_list = []   # [0, 1]  # individual sample name이 multiple vcf에서는 어디에 위치해 있는지
    for t in sample_name_individual_vcf:
        for s_i, s in enumerate (sample_name_multiple_vcf):    
            if t == s:
                s_list.append (s_i)
                break

    for line in ab:
        print ( "Multiple call - Individual call", line[0:5], sep = "\t")
        print ( *line[0:9], sep = "\t", end = "\t", file = output_file)
        for i in range ( len(s_list) - 1 ) :
            print (line [9 + s_list[i] ], sep = "\t", end = "\t", file = output_file)
        print (line [9 + s_list [ len(s_list) - 1 ] ], file = output_file)        
        
    output_file.close()
    output_file_unique.close()
    
    
