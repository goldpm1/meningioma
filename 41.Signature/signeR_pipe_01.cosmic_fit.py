import argparse
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity


parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--SIGNER_RESULT_Phat_PATH', type=str)
parser.add_argument('--SIGNER_RESULT_cosmic_fit_PATH', type=str)

args = parser.parse_args()


SIGNER_RESULT_Phat_PATH = args.SIGNER_RESULT_Phat_PATH
SIGNER_RESULT_cosmic_fit_PATH = args.SIGNER_RESULT_cosmic_fit_PATH 
COSMIC_PATH="/home/goldpm1/resources/COSMIC_SBS/COSMIC_v3.3.1_SBS_GRCh38.txt"

COSMIC_df = pd.read_csv (COSMIC_PATH, sep = "\t")
SIGNER_df = pd.read_csv (SIGNER_RESULT_Phat_PATH, sep = "\t")


# COSMIC_df의 signature 표현양식과 순서가 SIGNER_df와 다르니 맞춰주는 작업
new_COSMIC_df_index = []
for s in list (COSMIC_df ["Type"] ):
    t = s[0] + s[2] + s[-1] + ">" + s[4]
    new_COSMIC_df_index.append(t)

COSMIC_df ["Type"] = pd.Series(new_COSMIC_df_index)

# Create a dictionary to map the old row labels to the new row labels
row_mapping = {s: i for i, s in enumerate( list ( SIGNER_df.index ) )}

# Apply the mapping and reindex the DataFrame
COSMIC_df = COSMIC_df.set_index("Type").rename(index=row_mapping).sort_index().reset_index()

# 다시 index로 넣어주고 Type column은 제거해주기
COSMIC_df = COSMIC_df.set_index( SIGNER_df.index ).drop ("Type", axis = 1)


# Define two matrices
matrix1 = np.array([[1, 2, 3, 4 ], [3, 4, 6, 10], [5, 6, 10, 21]])    # SBS 개수 * 96 signatrues
matrix2 = np.array( [ [1, 2.2, 3.1, 4.4] ] )   # 1 * 96 signatures

# Compute cosine similarity
cos_sim = cosine_similarity( COSMIC_df.values.transpose(), SIGNER_df.values.transpose() )

# 가장 비슷한 (= cosine similarity가 가장 높은 SBS) SBS 찾기

output_file = open (SIGNER_RESULT_cosmic_fit_PATH, "w")
for signature_i in range ( cos_sim.shape[1] ):
    lst = sorted(  cos_sim[:, signature_i], reverse=True)
    lst_arg = []
    for i in lst:
        lst_arg.append ( list(cos_sim[:, signature_i]).index (i) )

    max = np.max ( cos_sim[:, signature_i] )
    max_arg = np.argmax ( cos_sim[:, signature_i] )
    #print ("S{}\tmax = {}\tmax_arg = {}".format(signature_i + 1, round (max, 2), COSMIC_df.columns[ max_arg ] )) 
    print ("S{}\tlst[0:5] = {}\tlst_arg[0:5] = {}".format (signature_i + 1, np.round (lst[0:5], 2) ,  COSMIC_df.columns [ [i for i in  lst_arg [0:5]] ]  ) , file = output_file )

output_file.close()
