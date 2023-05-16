import argparse
import gzip

parser = argparse.ArgumentParser(description='Here is usage direction.')
parser.add_argument('--INPUT_VCF_GZ', default="")
parser.add_argument('--OUTPUT_VCF', default="")

args = parser.parse_args()

INPUT_VCF_GZ = args.INPUT_VCF_GZ
OUTPUT_VCF = args.OUTPUT_VCF


input_file = gzip.open(INPUT_VCF_GZ, "r")
output_file = open(OUTPUT_VCF, "w")


chr_list = []
for idx in range(1, 23):
    chr_list.append("chr"+str(idx))
    chr_list.append(str(idx))
chr_list.append("chrX")
chr_list.append("chrY")


check = 0
while True:
    line = input_file.readline().decode('utf-8').rstrip('\n')

    if len(line) < 2:
        break

    if check == 0:   # Header를 출력
        print(line, file = output_file)

    line = line.split("\t")

    if line[0] in chr_list:
        check = 1
        print("\t".join(line), file=output_file)

input_file.close()
output_file.close()

