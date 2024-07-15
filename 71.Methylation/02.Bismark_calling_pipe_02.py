import os, subprocess, argparse
import pandas as pd

kwargs = {}

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument('--CALL_PATH', type=str)
parser.add_argument('--DSS_PATH', type=str)

args = parser.parse_args()

kwargs["CALL_PATH"] = args.CALL_PATH
kwargs["DSS_PATH"] = args.DSS_PATH

print (kwargs["CALL_PATH"])

df = pd.read_csv ( kwargs["CALL_PATH"], sep = "\t", header = None)
df.columns = [ "chr", "pos", "pos2", "methyl%", "X", "Non-methyl"]
df['N'] = df["X"] + df["Non-methyl"]

df = df.loc [ :, ["chr", "pos", "N", "X"]]
df = df[~df['chr'].str.contains('Un|decoy|alt|random|HLA')].reset_index (drop = True)

print (df)

df.to_csv ( kwargs["DSS_PATH"], sep = "\t", index = False)