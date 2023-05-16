import pandas as pd
import argparse

parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--Sample_ID', type=str, default="220930")
parser.add_argument('--TISSUE', type=str, default="Tumor")
parser.add_argument('--SEQUENZA_TO_PYCLONEVI_MATRIX_PATH', type=str, default="")
parser.add_argument('--FACETCNV_TO_PYCLONEVI_MATRIX_PATH', type=str, default="")

args = parser.parse_args()

Sample_ID = args.Sample_ID
TISSUE = args.TISSUE
SEQUENZA_TO_PYCLONEVI_MATRIX_PATH = args.SEQUENZA_TO_PYCLONEVI_MATRIX_PATH
FACETCNV_TO_PYCLONEVI_MATRIX_PATH = args.FACETCNV_TO_PYCLONEVI_MATRIX_PATH



############SEQUENZA RESULTS#############

# 읽고 chr_pos 로 정렬해줌
df = pd.read_csv (SEQUENZA_TO_PYCLONEVI_MATRIX_PATH, sep = "\t", names =["mutation_id", "sample_id", 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn', 'minor_cn', 'tumour_content', 'types' , 'gene', 'variant_classification'] )
#df = df.drop ('types', axis = 1)
df = df.drop_duplicates(["mutation_id", "sample_id"],  keep = 'first')

df.sort_values ( ['mutation_id', 'sample_id'], axis = 0, ascending = True, inplace = True)
df.to_csv (SEQUENZA_TO_PYCLONEVI_MATRIX_PATH, sep = "\t", index = False)




############FACETCNV RESULTS#############


df = pd.read_csv (FACETCNV_TO_PYCLONEVI_MATRIX_PATH, sep = "\t", names =["mutation_id", "sample_id", 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn', 'minor_cn', 'tumour_content', 'types', 'gene', 'variant_classification'] )
df = df.drop_duplicates(["mutation_id", "sample_id"],  keep = 'first')

df.sort_values ( ['mutation_id', 'sample_id'], axis = 0, ascending = True, inplace = True)
df.to_csv (FACETCNV_TO_PYCLONEVI_MATRIX_PATH, sep = "\t", index = False)
