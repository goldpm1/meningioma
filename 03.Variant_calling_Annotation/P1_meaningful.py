import pandas as pd
import numpy as np
import math
import os, vcf
import argparse

parser = argparse.ArgumentParser( description='The below is usage direction.')
parser.add_argument('--INPUT_VCF', type=str, default="")
parser.add_argument('--OUTPUT_DF', type=str, default="")
parser.add_argument('--AM_PATHOGENICITY', type=str, default="")
parser.add_argument('--SPLICEAI_SCORE', type=str, default="")

args = parser.parse_args()

vcf_reader = vcf.Reader(open( args.INPUT_VCF, "r"))
CSQ_title = vcf_reader.infos["CSQ"].desc.split ("Format: ")[1].split ("|")        # CSQ의 이름 
samplenames = vcf_reader.samples
samplenames_dict = { i : samplenames [i] for i in range ( len (samplenames ))  }
samplenames_dict_rev = { samplenames [i] : i for i in range ( len (samplenames ))  }

print (samplenames_dict)
df = pd.DataFrame ( columns = samplenames + ["Consequence", "am_pathogenicity", "SpliceAI_score", "Gene_Symbol", "PAF_max" ])

def POPULATION ( VEP_dict, info_dict, **kwargs ):
    for pop in kwargs ["population_db"]:
        if pop in VEP_dict.keys():
            if VEP_dict [ pop ] not in  ["", "."]:
                info_dict [ pop ] = float ( VEP_dict[ pop ] )
                if info_dict [ pop ]  > info_dict [ "PAF_max" ]:
                    info_dict [ "PAF_max" ] = info_dict [ pop ]
    return info_dict

def SPLICEAI ( VEP_dict ):
    filtered_list = [ float(value) for value in [VEP_dict["SpliceAI_pred_DS_AG"], VEP_dict["SpliceAI_pred_DS_AL"],  VEP_dict["SpliceAI_pred_DS_DG"], VEP_dict["SpliceAI_pred_DS_DL"]] if value != '' and value != None]
    if len ( filtered_list ) != 0:
        SpliceAI_score = np.max ( np.array ( filtered_list ) )
    else:
        SpliceAI_score = math.nan
    return SpliceAI_score

def ADD_MATRIX ( record, add):
    add_matrix = []
    for samplename in samplenames:
        add_matrix.append ( record.samples [ samplenames_dict_rev[samplename] ].data.AD )
    add_matrix = add_matrix + add
    return add_matrix

def main ( **kwargs ):
    line_num = 0
    for record in vcf_reader:        # record.CHROM, recrod.POS ,record.ALT
        CHR, POS, REF, ALT = record.CHROM, 	record.POS, record.REF,  record.ALT

        # if record.samples [ samplenames_dict_rev [ "240903_Blood" ] ].data.AD[1] != 0:    # Blood에서 alt >= 1이면 넘어가자
        #     continue

        info_dict = {}
        for pop in kwargs ["population_db"]:
            info_dict [ pop ] = "."
        info_dict [ "PAF_max"] = 0

        # df에 덧붙여줄 것
        global add_df
        add_df = pd.DataFrame ( columns = df.columns )


        u_list = []
        for u_index, u in enumerate( record.INFO["CSQ"] ) :   # Transcript 마다 돌기
            VEP_dict = {}
            for v_index, v in enumerate( u.split ("|") ) :
                VEP_dict [ CSQ_title [v_index] ] = v

            # Gene 정보
            GENE = VEP_dict ["SYMBOL"]

            # Population 정보
            info_dict = POPULATION ( VEP_dict, info_dict, **kwargs )

            # Splice_AI
            SpliceAI_score = SPLICEAI ( VEP_dict )
            if SpliceAI_score != math.nan:
                if SpliceAI_score > kwargs ["SpliceAI_score"]:
                    VEP_dict [ "Consequence" ] = "Splice_Site"
                    u_list.append ( u_index )
                    add_matrix = ADD_MATRIX ( record, [ VEP_dict ["Consequence"], math.nan, SpliceAI_score, GENE, info_dict ["PAF_max"] ]  )
                    add_df.loc [ len(add_df) ] = add_matrix
                    #print ( "\t{} : {}\t{}\t{}\t{}".format ( u_index, VEP_dict ["BIOTYPE"], VEP_dict ["Consequence"], VEP_dict ["am_pathogenicity"], SpliceAI_score  ) )


            # Canonical coding region
            if ( VEP_dict ["BIOTYPE"] == "protein_coding" ) & ( VEP_dict ["Consequence"] in  "missense_variant" ):
                if ( VEP_dict ["am_pathogenicity"] != "" ):
                    if ( float ( VEP_dict ["am_pathogenicity"] ) > kwargs ["am_pathogenicity"] ) :
                        u_list.append ( u_index )
                        add_matrix = ADD_MATRIX ( record, [ VEP_dict ["Consequence"], float ( VEP_dict ["am_pathogenicity"] ), math.nan, GENE, info_dict ["PAF_max"] ]  )
                        add_df.loc [ len(add_df) ] = add_matrix
                        #print ( "\t{} : {}\t{}\t{}\t{}".format ( u_index, VEP_dict ["BIOTYPE"], VEP_dict ["Consequence"], VEP_dict ["am_pathogenicity"], VEP_dict [ "SpliceAI_score" ]  ) )

            if ( VEP_dict ["BIOTYPE"] == "protein_coding" ) & ( VEP_dict ["Consequence"] in  kwargs ["Consequence"] ):
                u_list.append ( u_index )
                add_matrix = ADD_MATRIX ( record, [ VEP_dict ["Consequence"], math.nan, math.nan, GENE, info_dict ["PAF_max"] ]  )
                add_df.loc [ len(add_df) ] = add_matrix



        if u_list != []:
            add_df ["Consequence"] = pd.Categorical( add_df["Consequence"], categories = ["Splice_Site"] + kwargs ["Consequence"] + ["missense_variant"], ordered=True )
            add_df.sort_values(by = ["Consequence", "am_pathogenicity"] )

            # Choose most impactful transcript
            df.loc[ CHR + ":" + str(POS) ] = add_df.iloc [0, ]

            
            line_num += 1

        # if line_num >= 10:
        #     break
    print ( line_num )



kwargs = { "am_pathogenicity" : float(args.AM_PATHOGENICITY), "Consequence" : ["stop_gained", "frameshift_variant"], "SpliceAI_score" : float(args.SPLICEAI_SCORE), "population_db" : ["dbSNP_K1", "dbSNP_KRG", "dbSNP_EAS_AF", "gnomAD_AF" ] }

main ( **kwargs )

print ( df )

df.to_csv ( args.OUTPUT_DF, sep = "\t", index = True )