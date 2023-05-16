def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))


def visualization_decomposition( df, OUTPUT_SUPTITLE, ax ):
    ax.set_title(OUTPUT_SUPTITLE, fontsize = 24, fontweight='bold')
    ax.text(0.5, 0.9, "Purity of Tumor = {}\nPurity of Dura = {}".format(df[df['sample_id'].str.contains('Tumor')].iloc[0]["tumour_content"], df[df['sample_id'].str.contains('Dura')].iloc[0]["tumour_content"] ), ha='center', fontsize = 14 )
    ax.set_xlabel("VAF_Dura", fontdict = {"fontsize" : 14})
    ax.set_ylabel("VAF_Tumor", fontdict = {"fontsize" : 14})
    plt.style.use("seaborn-white")

    ax.set_xlim([0,  1])
    ax.set_ylim([0,  1])

    tabl = palettable.tableau.Tableau_20.mpl_colors
    vivid_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    colorlist = [i for i in vivid_10]
    shapelist = ["o", "s", "^", "v", ">", "<"]


    k = 0
    set_cluster_id = set([])

    while k < df.shape[0] - 1:
        i = 0 if df.iloc[k]["types"] == "Mutect2" else 3
        if df.iloc[k]["mutation_id"] == df.iloc[k + 1]["mutation_id"]:       # 2-shared mutation인 경우
            vaf_Dura = round ( int (df.iloc[k]["alt_counts"]) / ( int (df.iloc[k]["ref_counts"]) +  int (df.iloc[k]["alt_counts"])) , 2 )
            vaf_Tumor = round ( int (df.iloc[k + 1]["alt_counts"]) / ( int (df.iloc[k + 1]["ref_counts"]) +  int (df.iloc[k + 1 ]["alt_counts"])) , 2 )
            if df.iloc[k]["types"] == "Mutect2":
                print (vaf_Dura, vaf_Tumor, df.iloc[k]["mutation_id"], df.iloc[k]["gene"], df.iloc[k]["variant_classification"], sep = "\t")
            k2 = k + 2
        else:
            if "Dura" in df.iloc[k]["sample_id"]:
                vaf_Tumor = 0
                vaf_Dura = round ( int (df.iloc[k]["alt_counts"]) / ( int (df.iloc[k]["ref_counts"]) +  int (df.iloc[k]["alt_counts"])) , 2 )
            elif "Tumor" in df.iloc[k]["sample_id"]:
                vaf_Tumor = round ( int (df.iloc[k]["alt_counts"]) / ( int (df.iloc[k]["ref_counts"]) +  int (df.iloc[k]["alt_counts"])) , 2 )
                vaf_Dura = 0
            k2 = k + 1
            
        ax.scatter (vaf_Dura, vaf_Tumor,  alpha = 0.7, color = colorlist [ df.iloc[k]["cluster_id"]], marker = shapelist[i], s = 100, label = df.iloc[k]["cluster_id"] if df.iloc[k]["cluster_id"] not in set_cluster_id else ""   )
        set_cluster_id.add ( df.iloc[k]["cluster_id"]  )

        #legend_without_duplicate_labels ( ax )
        k = k2
    
    ax.legend()
    print ("\n\n")
        
        

if __name__ == "__main__":
    import pandas as pd
    import matplotlib.pyplot as plt
    import palettable, argparse
    import numpy as np


    parser = argparse.ArgumentParser( description='The below is usage direction.')
    parser.add_argument('--Sample_ID', type=str, default="220930")
    parser.add_argument('--SEQUENZA_TO_PYCLONEVI_MATRIX_PATH', type=str, default="")
    parser.add_argument('--SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH', type=str, default="")
    parser.add_argument('--FACETCNV_TO_PYCLONEVI_MATRIX_PATH', type=str, default="")
    parser.add_argument('--FACETCNV_TO_PYCLONEVI_OUTPUT_PATH', type=str, default="")
    parser.add_argument('--OUTPUT_PATH_SHARED', type=str, default="")
    parser.add_argument('--OUTPUT_PATH_TOTAL', type=str, default="")

    args = parser.parse_args()

    Sample_ID = args.Sample_ID
    SEQUENZA_TO_PYCLONEVI_MATRIX_PATH = args.SEQUENZA_TO_PYCLONEVI_MATRIX_PATH
    SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH = args.SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH
    FACETCNV_TO_PYCLONEVI_MATRIX_PATH = args.FACETCNV_TO_PYCLONEVI_MATRIX_PATH
    FACETCNV_TO_PYCLONEVI_OUTPUT_PATH = args.FACETCNV_TO_PYCLONEVI_OUTPUT_PATH
    OUTPUT_PATH_SHARED = args.OUTPUT_PATH_SHARED
    OUTPUT_PATH_TOTAL = args.OUTPUT_PATH_TOTAL


    df_seq_to_pycl = pd.read_csv (SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH, sep = "\t")
    df_seq_to_pycl = df_seq_to_pycl.drop_duplicates (['mutation_id'], keep = 'first').sort_values ( ['mutation_id'], axis = 0, ascending = True)[ [ "mutation_id", "cluster_id"] ].reset_index().drop ('index', axis = 1)
    df_facet_to_pycl = pd.read_csv (FACETCNV_TO_PYCLONEVI_OUTPUT_PATH, sep = "\t")
    df_facet_to_pycl = df_facet_to_pycl.drop_duplicates (['mutation_id'], keep = 'first').sort_values ( ['mutation_id'], axis = 0, ascending = True)[ [ "mutation_id", "cluster_id"] ].reset_index().drop ('index', axis = 1)

    df_seq_matrix = pd.read_csv (SEQUENZA_TO_PYCLONEVI_MATRIX_PATH, sep = "\t")
    df_facet_matrix = pd.read_csv (FACETCNV_TO_PYCLONEVI_MATRIX_PATH, sep = "\t")


    fig, ax = plt.subplots( figsize=(13, 6), nrows = 1, ncols = 2 )
    visualization_decomposition ( pd.concat ( [ df_seq_matrix, pd.Series ([1] * df_seq_matrix.shape[0] , name = "cluster_id")], axis = 1), "{} - Sequenza".format(Sample_ID), ax[0] )
    visualization_decomposition ( pd.concat ( [ df_facet_matrix, pd.Series ([1] * df_facet_matrix.shape[0] , name = "cluster_id")], axis = 1), "{} - facetcnv".format(Sample_ID), ax[1] )
    plt.savefig ( OUTPUT_PATH_TOTAL )
    plt.show()
    
    
    fig, ax = plt.subplots( figsize=(13, 6), nrows = 1, ncols = 2 )
    df_seq_integrated = pd.merge (df_seq_matrix, df_seq_to_pycl, left_on = "mutation_id", right_on = "mutation_id")
    visualization_decomposition ( df_seq_integrated, "{} - Sequenza".format(Sample_ID), ax[0] )
    df_facet_integrated = pd.merge (df_facet_matrix, df_facet_to_pycl, left_on = "mutation_id", right_on = "mutation_id")
    visualization_decomposition ( df_facet_integrated, "{} - facetcnv".format(Sample_ID), ax[1] )
    plt.savefig ( OUTPUT_PATH_SHARED )
    plt.show()