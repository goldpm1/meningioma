import pandas as pd
import numpy as np
import csv
import gzip
import glob
import time
import sys

MForiginal = sys.argv[1]

pd.set_option('display.max_seq_items', None)
pd.set_option('display.max_columns', None)

MForiginal_df = pd.read_csv(MForiginal, delim_whitespace=True)


MForiginal_df["chr"] = MForiginal_df["id"].str.split("~").str[1]
MForiginal_df["start_pos"] = MForiginal_df["id"].str.split("~").str[2].astype(int); MForiginal_df["start_pos"] = MForiginal_df["start_pos"] - 1
MForiginal_df["end_pos"] = MForiginal_df["id"].str.split("~").str[2]
MForiginal_df["ref"] = MForiginal_df["id"].str.split("~").str[3]
MForiginal_df["alt"] = MForiginal_df["id"].str.split("~").str[4]

MForiginal_df[["chr","start_pos","end_pos","ref","alt"]].to_csv (MForiginal + ".bed", index = False, na_rep = '.', sep = '\t', header = False)
