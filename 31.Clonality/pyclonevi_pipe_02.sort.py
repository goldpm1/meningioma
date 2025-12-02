import pandas as pd
import numpy as np
import argparse
import warnings
warnings.simplefilter (action = 'ignore')

parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--INPUT_TSV', type=str, default="/data/project/Meningioma/31.Clonality/02.pyclonevi/230303/230303.facetcnv_to_pyclonevi.tsv")
parser.add_argument('--OUTPUT_TSV', type=str, default="/data/project/Meningioma/31.Clonality/02.pyclonevi/230303/230303.facetcnv_to_pyclonevi.tsv")

args = parser.parse_args()

kwargs = {}
INPUT_TSV=args.INPUT_TSV
OUTPUT_TSV=args.OUTPUT_TSV



df = pd.read_csv ( INPUT_TSV, sep = "\t")
df["axis"] = df["sample_id"].str.split ("_").str[-1]

df2 = pd.DataFrame ( df.groupby( ["cluster_id", "axis"])["cellular_prevalence"].mean() ).unstack(fill_value=0).reset_index()

samplename_list = []
for i in range (1, len(df2.columns) ):
    samplename_list.append ( df2.columns[i][1] )
df2.columns = ["cluster_id"] + samplename_list

df2 = df2.sort_values(by = ["Tumor", "Dura"], ascending=False).reset_index()
df2["index"] = df2.index

# 새로운 index를 부여하기 위해서 merge함
new_df = pd.merge (df, df2, left_on  = "cluster_id", right_on= "cluster_id", how='left')
selected_columns = [col for col in new_df.columns if col not in ["cluster_id"] + samplename_list ]
new_df = new_df [selected_columns]

# 새로운 index를 cluster_df로 바꿔줌
new_df.rename(columns = {"index":"cluster_id" }, inplace = True)
new_df = new_df.sort_values(by = ["cluster_id", "mutation_id"], ascending=True).reset_index ( drop = True)

new_df.to_csv (OUTPUT_TSV, sep = "\t", index = False)
print (new_df)