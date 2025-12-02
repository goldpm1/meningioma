import argparse
import math
import numpy as np

# Constants
DEFAULT_DIVISOR = 14906 ##samples in ./GENIE/data_NSCLC_all.txt (GENIE v15)
ASIAN_DIVISOR = 1069 ##samples in ./GENIE/data_NSCLC_asian.txt (GENIE v15)
NON_ASIAN_DIVISOR = 13837 ##samples in ./GENIE/data_NSCLC_nonasian.txt (GENIE v15)
DEFAULT_MUTATION_FREQUENCY = 10**-6
IPM_PRIOR = 0.01  ##사전정보가 없어서 0.5
MPLC_PRIOR = 0.99 ##사전정보가 없어서 0.5

# Function to read mutation frequencies
#def read_mutation_frequencies(filename, divisor):
#   mutation_frequencies = {}
#    with open(filename, 'r') as file:
#        for line in file:
#            gene, hgvsc, count = line.split('\t')
#            ID = f'{gene}\t{hgvsc}'
#            mutation_frequencies[ID] = float(count) / divisor
#    return mutation_frequencies

# Function to parse a line of input file
def parse_line(line):
    fields = line.rstrip().split('\t')
    gene, hgvsp, avaf, bvaf = fields[1], fields[3], float(fields[4]), float(fields[5])
    return gene, hgvsp, avaf, bvaf

# Function to calculate the probability
def calculate_probability(avaf, bvaf, mutation_frequency):
    if avaf >= 0.5:
        avaf = 0.49999
    if  bvaf >= 0.5:
        bvaf = 0.49999

    if avaf != 0 and bvaf != 0:
        return math.log10((2 * avaf + (1 - 2 * avaf) * mutation_frequency) / mutation_frequency)
    elif avaf != 0 and bvaf == 0:
        return math.log10(1 - 2 * (avaf))
    else:
        return 0  # Log(1) is 0

# Function to calculate score and confidence level
def calculate_score_and_confidence(prob_list, mutation_prior_ratio, max_mutation_count):
    score = sum(prob_list) + math.log10(mutation_prior_ratio)
    classification = 'IPM' if score > 0 else 'MPLC'
    confidence_level = get_confidence_level(max_mutation_count)
    return score, classification, confidence_level

# Function to determine confidence level based on maximum mutation count
def get_confidence_level(max_mutation_count):
    if max_mutation_count <= 2:
        return 'Likely'
    else:
        return 'Confident'


# Main processing function
def process_metel(input_file, mode, race):
    with open(input_file, 'r') as file:
        #next(file)  # Skip the header line
        ab_list, ba_list = [], []
        a_mutation_count, b_mutation_count = 0, 0  

        for line in file:
            gene, hgvsp, avaf, bvaf = parse_line(line)
            if gene == 'NF2' :
               mutation_frequency = 0.8  ####NF2의 mutation frequency
               ab_list.append(calculate_probability(avaf , bvaf, mutation_frequency))
               if mode == 'syn':
                   ba_list.append(calculate_probability(bvaf , avaf, mutation_frequency))
                   print ( "{}\tmutation frequency = {}\tavaf = {}\tbvaf = {}\tab_list = {}\tba_list = {}".format ( gene, mutation_frequency, round ( avaf / 2, 3 ), round ( bvaf, 3 ), np.round ( calculate_probability(avaf / 2, bvaf, mutation_frequency), 3 ), np.round ( calculate_probability(bvaf / 2, avaf, mutation_frequency), 3 ) ) )
            else :
                mutation_frequency = DEFAULT_MUTATION_FREQUENCY
                if (avaf > 0)  & (bvaf > 0) & (avaf < 0.2):       # subclonal shared mutation의 경우
                    mutation_frequency = 0.6
                if  gene in ["TRAF7", "KLF4", "AKT1"]:
                    mutation_frequency = 0.4
                ab_list.append(calculate_probability(avaf, bvaf, mutation_frequency))
                if mode == 'syn':
                    ba_list.append(calculate_probability(bvaf, avaf, mutation_frequency))
                    print ( "{}\tmutation frequency = {}\tavaf = {}\tbvaf = {}\tab_list = {}\tba_list = {}".format ( gene, mutation_frequency, round ( avaf, 3 ), round ( bvaf, 3 ), np.round ( calculate_probability(avaf, bvaf, mutation_frequency), 3 ), np.round ( calculate_probability(bvaf, avaf, mutation_frequency), 3 ) ) )
            
            # Update mutation counts for A and B
            if avaf != 0:
                a_mutation_count += 1
            if bvaf != 0:
                b_mutation_count += 1


        print ( "\n\tab_list  (IPM, len = {} ) = {}".format ( len ( ab_list ), np.round ( ab_list, 3 ) ) )
        print ( "\tba_list  (MPLC, len = {} ) = {}".format ( len ( ba_list ), np.round ( ba_list, 3 ) ) )

        max_mutation_count = max(a_mutation_count, b_mutation_count)
        mutation_prior_ratio = IPM_PRIOR / MPLC_PRIOR

        score_ab, classification_ab, confidence_ab = calculate_score_and_confidence(ab_list, mutation_prior_ratio, max_mutation_count)
        score_ba, classification_ba, confidence_ba = calculate_score_and_confidence(ba_list, mutation_prior_ratio, max_mutation_count)

        print ( "\n\tscore_ab = {}\tscore_ba = {}".format ( np.round ( score_ab, 3 ), np.round ( score_ba, 3 ) ) )

        # Selecting the result with the higher absolute score
        if abs(score_ab) > abs(score_ba):
            final_score, final_classification, final_confidence = score_ab, classification_ab, confidence_ab
        else:
            final_score, final_classification, final_confidence = score_ba, classification_ba, confidence_ba
        
        print ( "\n\tfinal_score = {}\tfinal_classification = {}\tfinal_confidence = {}".format ( round ( final_score, 3 ), final_classification, final_confidence) )
        #print ( "Classification_Score(s)\tDiagnosis_Result\tConfidence_Level\tRace" )


        return final_score, final_classification, final_confidence, race






#################################################################



# Setup argument parser
parser = argparse.ArgumentParser(description='Process VAF data for MeTel algorithm.')
parser.add_argument('input_file', type=str, help='Input file name')
parser.add_argument('-s', '--synmeta', type=str, default='syn', choices=['syn', 'meta'], help='Synchronicity mode')
parser.add_argument('-r', '--race', type=str, choices=['asian', 'non-asian'], help='Race mode')
parser.add_argument('output_file', type=str, help='Output file name')

# Parse arguments
args = parser.parse_args()

# Determine mutation frequency file and divisor based on race option
#mutation_freq_file = './GENIE/data_mutations_c.txt'  # Default file
#divisor = DEFAULT_DIVISOR

#if args.race == 'asian':
#    mutation_freq_file = './GENIE/data_mutations_c_asian.txt'
#    divisor = ASIAN_DIVISOR
#elif args.race == 'non-asian':
#    mutation_freq_file = './GENIE/data_mutations_c_nonasian.txt'
#    divisor = NON_ASIAN_DIVISOR

# Read mutation frequencies
#mutation_frequencies = read_mutation_frequencies(mutation_freq_file, divisor)
 

# Determine race value for output
race_value = 'Unspecified'
if args.race:
    race_value = args.race

# Process input file
result = process_metel(args.input_file, args.synmeta, race_value)

# Write the result to output file
with open(args.output_file, 'w') as out:
    out.write('Classification_Score(s)\tDiagnosis_Result\tConfidence_Level\tRace\n')
    out.write('\t'.join(map(str, result)) + '\n')
