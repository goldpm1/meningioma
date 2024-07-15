import os, subprocess, argparse
import pandas as pd
import numpy as np
import vcf, pysam
from scipy.stats import binom

kwargs = {}

parser = argparse.ArgumentParser(description='The below is usage direction.')
# parser.add_argument('--BCFTOOLS_MERGE_OUTPUT_VCF', type=str, default = "/data/project/Meningioma/61.Lowinput/01.XT_HS/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/190426.merge.vcf")
# parser.add_argument('--OUTPUT_HEATMAP_PATH', type=str, default = "/data/project/Meningioma/61.Lowinput/01.XT_HS/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/190426.heatmap.pdf" )

parser.add_argument('--BCFTOOLS_MERGE_OUTPUT_VCF', type=str, default = "/data/project/Meningioma/61.Lowinput/02.PTA/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/230405.merge.vcf")
parser.add_argument('--OUTPUT_HEATMAP_PATH', type=str, default = "/data/project/Meningioma/61.Lowinput/02.PTA/06.HC/07.2D_merged/01.BCFTOOLS_MERGE_TXT/230405.heatmap.pdf" )

args = parser.parse_args()

kwargs["BCFTOOLS_MERGE_OUTPUT_VCF"] = args.BCFTOOLS_MERGE_OUTPUT_VCF 
kwargs["OUTPUT_HEATMAP_PATH"] = args.OUTPUT_HEATMAP_PATH

print ( kwargs )


####################################################################################################################################


def find_not_blood_sample(samplenames):
  li = []
  for i, name in enumerate(samplenames):
    if "Blood" not in name:
        li.append (i)
  return li


import pandas as pd
import vcf

vcf_reader = vcf.Reader(open( kwargs["BCFTOOLS_MERGE_OUTPUT_VCF"], "r"))
samplenames = vcf_reader.samples
not_blood_sample_i = find_not_blood_sample ( samplenames )

df = pd.DataFrame ( columns = [ samplenames[i] for i in not_blood_sample_i ] )            # Blood가 아닌 샘플들에 대해서만 집어넣겠다
chrpos = []

for line in vcf_reader:
    CHR = line.CHROM
    POS = str (line.POS)
    REF, ALT = str(line.REF), str( line.ALT[0] )

    matrix = []
    for sample_i, samplename in enumerate (samplenames):
        if line.samples [sample_i].data.GT in ["./.", "0/0", "0|0"]:
            matrix.append ( 0 )
        else:
            matrix.append ( 1  )

    df.loc [ len(df.index) ] =  [ matrix[i] for i in not_blood_sample_i ]      # 한줄씩 집어넣기
    chrpos.append ( CHR + "_" + POS )
df.index =  chrpos

df = df.sort_values(by=df.columns.tolist(), ascending=False)
print (df)
df.to_csv ( kwargs["OUTPUT_HEATMAP_PATH"].replace("pdf", "tsv"), sep = "\t")


############# Clone, subclone 별로 색깔 다륵 ㅔ해주기 ########################33
np_matrix = np.array (df)

for j in range ( np_matrix.shape[0] ) :
    if all ( np_matrix [j] == 1):
        np_matrix [j] = 1
    elif ( np_matrix[j][0] == 1 ) & ( all ( np_matrix [j][1:] == 0) ) :
        np_matrix [j][ 0 ] = 2
    elif ( np_matrix[j][0] == 0 ) & ( all ( np_matrix [j][1:] == 1) ) :
        np_matrix [j][ 1: ] = 3
    else:
        for i in range (np_matrix.shape[1]):
            if np_matrix [j] [i]  == 1:
                np_matrix [j][i] = 4




#def HEATMAP_VISUALIZATION (df, title, Output_filename, **kwargs):
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcl
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap

title = kwargs["OUTPUT_HEATMAP_PATH"].split("/")[-1]
Output_filename = kwargs["OUTPUT_HEATMAP_PATH"]

plt.rcParams["font.family"] = 'arial'

# Define the colors
fig, ax = plt.subplots ( nrows = 1, ncols = 1, figsize =(25 , 6 ))

# Create the colormap
cmap_dict = { 0: "white", 1: "#D7FFA4", 2: "#477045", 3: "#4D63DD", 4: "#5AA0FF" }

fig.subplots_adjust ( wspace = 0.4, bottom = 0.03, top = 0.7, left = 0.22, right = 0.98)
fig.set_facecolor('white')
sns.heatmap ( np_matrix.T, cmap = ListedColormap([ cmap_dict[i] for i in range(5)])  , linewidths = 0, linecolor = "black" )   # fmt=".2f", 
    
    
fig.suptitle ( title, fontsize = 12, fontweight = "bold", ha = "left", x = 0.5 )
ax.set_xticklabels( [ i for i in df.index ] , fontsize = 12, ha = 'left' )
ax.tick_params(axis = 'x',  rotation = 45)
ax.set_yticklabels( [ i.split("_")[-1] for i in df.columns], fontsize = 12, va = 'center' )
ax.tick_params(axis = 'y',  rotation = 0)


plt.tick_params(axis = 'both', which = 'major', labelsize = 12, left = False, labelbottom = False, bottom=False, top = False, labeltop=True)

fig.savefig ( Output_filename, dpi = 300)