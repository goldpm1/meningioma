import argparse

parser = argparse.ArgumentParser(description='Here is usage direction.')
parser.add_argument('--INPUT', default="")
parser.add_argument('--OUTPUT', default="")

args = parser.parse_args()

INPUT = args.INPUT
OUTPUT = args.OUTPUT


input_file = open(INPUT, "r")
output_file = open(OUTPUT, "w")


chr_list = []
for idx in range(1, 23):
    chr_list.append("chr"+str(idx))
    chr_list.append(str(idx))
chr_list.append("chrX")
chr_list.append("chrY")
chr_list.append("chrM")


while True:
    line = input_file.readline().rstrip('\n')

    if len(line) < 2:
        break

    if line[0] == "#":
        print(line, file=output_file)
        continue

    line = line.split("\t")

    if line[0] in chr_list:
        print("\t".join(line), file=output_file)

input_file.close()
output_file.close()
