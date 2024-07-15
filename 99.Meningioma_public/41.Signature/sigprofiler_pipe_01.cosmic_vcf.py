import argparse
from SigProfilerAssignment import Analyzer as Analyze
import sigProfilerPlotting as sigPlt
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns


parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--Sample_ID', type=str)
parser.add_argument('--SIGPROFILER_INPUT_VCF_DIR', type=str)
parser.add_argument('--SIGPROFILER_RESULT_DIR', type=str)

args = parser.parse_args()

kwargs = {}
kwargs["Sample_ID"] = args.Sample_ID
kwargs["SIGPROFILER_INPUT_VCF_DIR"] = args.SIGPROFILER_INPUT_VCF_DIR
kwargs["SIGPROFILER_RESULT_DIR"] = args.SIGPROFILER_RESULT_DIR 


Analyze.cosmic_fit(kwargs["SIGPROFILER_INPUT_VCF_DIR"], kwargs["SIGPROFILER_RESULT_DIR"], input_type = "vcf", genome_build="GRCh38", 
                   context_type = "96", exome =True,
                   signatures = None, signature_database = None, collapse_to_SBS96=False, make_plots = True,
                   exclude_signature_subgroups = [ "Chemotherapy_signatures", "Treatment_signatures", "Artifact_signatures", "Test2"] )


# output file 읽기
f = open(kwargs["SIGPROFILER_RESULT_DIR"] + '/Assignment_Solution/Activities/Assignment_Solution_Activities.txt','r')

data = f.readlines()

column_name = data[0].rstrip().split('\t')
data_split = [x.rstrip().split('\t') for x in data[1:]]

df = pd.DataFrame(data_split, columns = column_name)
df = df.set_index(keys = 'Samples')

for i in df.columns:
    df[i] = pd.to_numeric(df[i])
    
    

# 0으로 가득찬 df 를 filter하기    
for column in df.columns:
    check = set()
    for value in df[column]:
        check.add(value)
    if len(check) == 1:
        df = df.drop(column, axis = 'columns')


# percentage로 변환하기
SBS_name = df.columns
total_count = {}

for index, row in df.iterrows():  ## index == A.~~ row = value
    total = 0
    for value in row:
        total += value
    total_count[index] = total

for index, row in df.iterrows():
    for idx, value in enumerate(row):
        df[SBS_name[idx]][index] = value/total_count[index]*100
        
df.to_csv (kwargs["SIGPROFILER_RESULT_DIR"] + "/results_df.tsv", sep = "\t")


# Stacked bar chart
fig, ax = plt.subplots( figsize=(12, 6), nrows = 1, ncols = 1 )
fig.subplots_adjust ( bottom = 0.25, top = 0.88, left = 0.04, right = 0.96)
sns.set_style("white")
ax.set_title ("Stacked bar chart", fontsize = 24, fontweight='heavy')
df.plot (kind = "bar", stacked=True , ax = ax)
ax.set_xticklabels (list ( df.index ), rotation  = 90, fontsize = 11, fontweight='semibold')
ax.set_xlabel("")
plt.savefig (kwargs["SIGPROFILER_RESULT_DIR"] + "/sigprofiler_barchart.jpg", dpi = 300)
plt.show()