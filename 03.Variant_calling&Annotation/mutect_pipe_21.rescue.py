import argparse
import numpy as np
import pandas as pd
import os, pysam, vcf, pybedtools
import gzip


parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--Sample_ID', type=str)
parser.add_argument('--TISSUE', type=str)
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
    output_file = open ( kwargs["INDIVIDUAL_RESCUED_VCF_PATH"] + ".unsorted" , "w" )
    output_file_unique = open ( kwargs["INDIVIDUAL_UNIQUE_VCF_PATH"]  , "w" )
    
    input_file = open ( kwargs["INDIVIDUAL_VCF_PATH"], "r")
    for line in input_file.readlines():
        line = line.rstrip("\n")            
        if line [0] == '#':        # Header 저장
            print (line, file = output_file)
            print (line, file = output_file_unique)
            continue
        else:
            break
    input_file.close()
    
    
    
    # Mutect call ∩ Individual call  → 그냥 Individual call을 출력
    ab = pybed_individual.intersect(pybed_multiple)
    for line in ab:
        print (*line, sep = "\t", file = output_file)
        
    # Individual call - Multiple call  → 그냥 Individual call을 출력
    ab = pybed_individual.intersect(pybed_multiple, v = True)
    for line in ab:
        print (*line, sep = "\t", file = output_file)
        print (*line, sep = "\t", file = output_file_unique)
        
        
    # Multiple call - Individual call → Multiple call을 이식해서 출력해줌
    ab = pybed_multiple.intersect(pybed_individual, v = True)
    s_list = []   # [0, 1]  # individual sample name이 multiple vcf에서는 어디에 위치해 있는지
    for t in sample_name_individual_vcf:
        for s_i, s in enumerate (sample_name_multiple_vcf):    
            if t == s:
                s_list.append (s_i)
                break

    for line in ab:
        print ( *line[0:9], sep = "\t", end = "\t", file = output_file)
        for i in range ( len(s_list) - 1 ) :
            print (line [9 + s_list[i] ], sep = "\t", end = "\t", file = output_file)
        print (line [9 + s_list [ len(s_list) - 1 ] ], file = output_file)        
        
    output_file.close()
    output_file_unique.close()
    
    
    print ( "/opt/Yonsei/bcftools/1.7/bcftools sort -Ov " + kwargs["INDIVIDUAL_RESCUED_VCF_PATH"] + ".unsorted" + " -o " + kwargs["INDIVIDUAL_RESCUED_VCF_PATH"] )
    os.system ("/opt/Yonsei/bcftools/1.7/bcftools sort -Ov " + kwargs["INDIVIDUAL_RESCUED_VCF_PATH"] + ".unsorted" + " -o " + kwargs["INDIVIDUAL_RESCUED_VCF_PATH"])
    os.system ("rm -rf " + kwargs["INDIVIDUAL_RESCUED_VCF_PATH"] + ".unsorted")    