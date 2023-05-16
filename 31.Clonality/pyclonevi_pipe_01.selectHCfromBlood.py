import pandas as pd
import random
import argparse
import pybedtools

parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--RANDOM_PICK', type=int, default=50)
parser.add_argument('--HC_OUTPUT_PATH', type=str, default="/home/goldpm1/Meningioma/06.hc/03.HF/220930/Blood/220930_Blood.DP100.vcf")
parser.add_argument('--HC_BLOOD_RANDOM_PICK_PATH', type=str, default="/home/goldpm1/Meningioma/31.Clonality/01.make_matrix/220930/220930.random_pick_50.bed")

args = parser.parse_args()

kwargs = {}
kwargs["RANDOM_PICK"] = args.RANDOM_PICK
kwargs["HC_OUTPUT_PATH"] = args.HC_OUTPUT_PATH
kwargs["HC_BLOOD_RANDOM_PICK_PATH"] = args.HC_BLOOD_RANDOM_PICK_PATH

hc_pybed_object = pybedtools.BedTool( kwargs["HC_OUTPUT_PATH"] )
random_index = sorted ( random.sample (range (len( hc_pybed_object )), kwargs["RANDOM_PICK"]) )          # HC call 50개만 골라줌 → segment와 intersect할 경우 50개 미만이 될 것

output_file = open( kwargs["HC_BLOOD_RANDOM_PICK_PATH"], "w")
for k in random_index:
    interval = hc_pybed_object [k]
    #print (interval [0:5])
    CHR, POS = str(interval[0]), int(interval[1])
    print ("{}\t{}\t{}".format(CHR, str(POS - 1), str(POS + 1)), file = output_file)
output_file.close()