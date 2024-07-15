def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))


def visualization_decomposition( df, OUTPUT_SUPTITLE, ax ):
    ax.set_title(OUTPUT_SUPTITLE, fontsize = 11, fontweight='bold')
    #ax.text(0.5, 0.9, "Purity of Tumor = {}\nPurity of Dura = {}".format(df[df['sample_id'].str.contains('Tumor')].iloc[0]["tumour_content"], df[df['sample_id'].str.contains('Dura')].iloc[0]["tumour_content"] ), ha='center', fontsize = 14 )
    ax.set_xlabel("VAF_Dura", fontdict = {"fontsize" : 8}, labelpad = 2 )
    ax.set_ylabel("VAF_Tumor", fontdict = {"fontsize" : 8}, labelpad = 2 )
    ax.tick_params( axis = 'x', labelsize = 8, pad = -1 )
    ax.tick_params( axis = 'y', labelsize = 8, pad = -1 )
    sns.set_style("white") 
    for axis in ['left', 'right', 'top', 'bottom']:
        ax.spines[axis].set_linewidth ( 1.5 )




    vaf_Tumor_list = []
    vaf_Dura_list = []

    tabl = palettable.tableau.Tableau_20.mpl_colors
    vivid_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    if len ( set ( df["cluster_id"] ) ) >= 4:
        colorlist = ["#6AEF6C", "#BCA97B", "#8D493A", "#7288C6"] + [i for i in vivid_10]
    else:
        colorlist = ["#6AEF6C", "#BCA97B", "#7288C6"] + [i for i in vivid_10]
    shapelist = ["o", "s", "^", "v", ">", "<"]

    
    df[ ["CHR", "POS", "ALT", "REF"] ] = df["mutation_id"].str.split("_", expand = True)
    df["CHR"] = pd.Categorical(df["CHR"], 
                                                categories = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"], 
                                                ordered=True)
    df = df.sort_values (by = ['cluster_id', 'CHR', 'POS'], axis = 0).drop (["CHR", "POS"], axis = 1).reset_index(drop = True)

    df_count = np.unique ( df ["cluster_id"], return_counts = True )

    # for i in range ( len(df_count[0]) ) :
    #     ax.text( ax.get_xlim()[1] / 2, ax.get_xlim()[1] / 1.5 - 0.03* df_count[0][i], "cluster{} = {}".format( df_count[0][i], int(df_count[1][i] / 2)  ), ha = 'center', fontsize = 8 )

    k = 0
    set_cluster_id = set([])

    while k < df.shape[0] - 1:
        if df.iloc[k]["mutation_id"] == df.iloc[k + 1]["mutation_id"]:       # 2-shared mutation인 경우
            vaf_Dura = round ( int (df.iloc[k]["alt_counts"]) / ( int (df.iloc[k]["ref_counts"]) +  int (df.iloc[k]["alt_counts"])) , 3 )
            vaf_Tumor = round ( int (df.iloc[k + 1]["alt_counts"]) / ( int (df.iloc[k + 1]["ref_counts"]) +  int (df.iloc[k + 1 ]["alt_counts"])) , 3 )
            i = 3 if ( (int (df.iloc[k]["alt_counts"]) == 0) | ( int (df.iloc[k + 1]["alt_counts"]) == 0)) else 0
            k2 = k + 2
        else:
            i = 3
            if "Dura" in df.iloc[k]["sample_id"]:
                vaf_Tumor = 0
                vaf_Dura = round ( int (df.iloc[k]["alt_counts"]) / ( int (df.iloc[k]["ref_counts"]) +  int (df.iloc[k]["alt_counts"])) , 3 )
            elif "Tumor" in df.iloc[k]["sample_id"]:
                vaf_Tumor = round ( int (df.iloc[k]["alt_counts"]) / ( int (df.iloc[k]["ref_counts"]) +  int (df.iloc[k]["alt_counts"])) , 3 )
                vaf_Dura = 0
            k2 = k + 1

        size = [ 100, 150, 150, 150 ]
        ax.scatter ( vaf_Dura, vaf_Tumor,  alpha = 0.7, s = size[i], 
                        color = colorlist [ df.iloc[k]["cluster_id"]], 
                        marker = shapelist[i], linewidths=0,
                        label = "{}".format ( df.iloc [k]["cluster_id"]) )
        if df.iloc[k]["gene"] in ["NF2", "AKT1", "KLF4", "TRAF7"]:
            ax.text ( vaf_Dura, vaf_Tumor, df.iloc[k]["gene"],   ha = "left", va = "bottom", fontdict = {"fontsize": 11, "fontweight" : "bold", "fontstyle": "italic"} )
        vaf_Tumor_list.append (vaf_Tumor)
        vaf_Dura_list.append (vaf_Dura)
        set_cluster_id.add ( df.iloc[k]["cluster_id"]  )

        k = k2
    
    #legend_without_duplicate_labels ( ax )


    # ax 정하기

    if ("190426" in OUTPUT_SUPTITLE) | ("221102" in OUTPUT_SUPTITLE) | ("230127" in OUTPUT_SUPTITLE) | ("230405" in OUTPUT_SUPTITLE) | ("230419" in OUTPUT_SUPTITLE) | ("230822" in OUTPUT_SUPTITLE) | ("230920" in OUTPUT_SUPTITLE):
        ax.set_xlim([0,  0.7]);  ax.set_ylim([0,  0.7])
    elif "220930" in OUTPUT_SUPTITLE:
        ax.set_xlim([0,  0.55]);  ax.set_ylim([0,  0.55])
    else:
        ax.set_xlim([0,  0.45]); ax.set_ylim([0,  0.45])

    ax.text( ax.get_xlim()[1] / 2 , ax.get_ylim()[1] / 1.2, "Purity of Tumor = {}\nPurity of Dura = {}".format(df[df['sample_id'].str.contains('Tumor')].iloc[0]["tumour_content"], df[df['sample_id'].str.contains('Dura')].iloc[0]["tumour_content"] ), ha='center', fontsize = 8 )
    for i in range ( len(df_count[0]) ) :
        ax.text( ax.get_xlim()[1] / 2, ax.get_ylim()[1] / 1.2 - 0.03* (i + 1), "cluster{} = {}".format( df_count[0][i], int(df_count[1][i] / 2)  ), ha = 'center', fontsize = 8 )


    print ("\n\n")
    
    return df, ax
        





def visualization_decomposition_scaled ( df, OUTPUT_SUPTITLE, ax ):
    ax.set_title(OUTPUT_SUPTITLE, fontsize = 11, fontweight='bold')
    ax.set_xlabel("VAF_Dura", fontdict = {"fontsize" : 8}, labelpad = 2 )
    ax.set_ylabel("VAF_Tumor", fontdict = {"fontsize" : 8}, labelpad = 2 )
    ax.tick_params( axis = 'x', labelsize = 8, pad = -1 )
    ax.tick_params( axis = 'y', labelsize = 8, pad = -1 )
    sns.set_style("white") 
    for axis in ['left', 'right', 'top', 'bottom']:
        ax.spines[axis].set_linewidth ( 1.5 )


    if ("221102" in OUTPUT_SUPTITLE) | ("221202" in OUTPUT_SUPTITLE) | ("230303" in OUTPUT_SUPTITLE) | ("230419" in OUTPUT_SUPTITLE) | ("230526" in OUTPUT_SUPTITLE):
        MULTIPLIER = 5        # Dura 0.03 까지 보여준다
        xtick = np.round ( np.arange (0, 0.21, 0.01)  , 2) 
    elif ("230127" in OUTPUT_SUPTITLE) | ("230405" in OUTPUT_SUPTITLE) | ("230920" in OUTPUT_SUPTITLE)   :
        MULTIPLIER = 3       # Dura 0.05까지 보여준다
        xtick = np.round ( np.arange (0, 0.21, 0.01)  , 2)      
    elif ("230822" in OUTPUT_SUPTITLE)  :
        MULTIPLIER = 2       # 
        xtick = np.round ( np.arange (0, 0.21, 0.01)  , 2)      
    elif ("221026" in OUTPUT_SUPTITLE):
        MULTIPLIER = 0.75  # Dura 0.25까지 보여준다
        xtick = np.round ( np.arange (0, 0.3, 0.05)   , 2 )
    else:
        MULTIPLIER = 1    # Dura 0.18까지 보여준다
        xtick = np.round ( np.arange (0, 0.25, 0.05)  , 2)
    


    vaf_Tumor_list = []
    vaf_Dura_list = []

    tabl = palettable.tableau.Tableau_20.mpl_colors
    vivid_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    if len ( set ( df["cluster_id"] ) ) >= 4:
        colorlist = ["#6AEF6C", "#BCA97B", "#8D493A", "#7288C6"] + [i for i in vivid_10]
    else:
        colorlist = ["#6AEF6C", "#BCA97B", "#7288C6"] + [i for i in vivid_10]
    shapelist = ["o", "s", "^", "v", ">", "<"]

    
    df[ ["CHR", "POS", "ALT", "REF"] ] = df["mutation_id"].str.split("_", expand = True)
    df["CHR"] = pd.Categorical(df["CHR"], 
                                                categories = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"], 
                                                ordered=True)
    df = df.sort_values (by = ['cluster_id', 'CHR', 'POS'], axis = 0).drop (["CHR", "POS"], axis = 1).reset_index(drop = True)

    df_count = np.unique ( df ["cluster_id"], return_counts = True )

    k = 0
    set_cluster_id = set([])

    while k < df.shape[0] - 1:
        if df.iloc[k]["mutation_id"] == df.iloc[k + 1]["mutation_id"]:       # 2-shared mutation인 경우
            vaf_Dura = round ( int (df.iloc[k]["alt_counts"]) / ( int (df.iloc[k]["ref_counts"]) +  int (df.iloc[k]["alt_counts"])) , 3 )
            vaf_Tumor = round ( int (df.iloc[k + 1]["alt_counts"]) / ( int (df.iloc[k + 1]["ref_counts"]) +  int (df.iloc[k + 1 ]["alt_counts"])) , 3 )
            i = 3 if ( (int (df.iloc[k]["alt_counts"]) == 0) | ( int (df.iloc[k + 1]["alt_counts"]) == 0)) else 0
            k2 = k + 2
        else:
            i = 3
            if "Dura" in df.iloc[k]["sample_id"]:
                vaf_Tumor = 0
                vaf_Dura = round ( int (df.iloc[k]["alt_counts"]) / ( int (df.iloc[k]["ref_counts"]) +  int (df.iloc[k]["alt_counts"])) , 3 )
            elif "Tumor" in df.iloc[k]["sample_id"]:
                vaf_Tumor = round ( int (df.iloc[k]["alt_counts"]) / ( int (df.iloc[k]["ref_counts"]) +  int (df.iloc[k]["alt_counts"])) , 3 )
                vaf_Dura = 0
            k2 = k + 1
            
        print ( vaf_Dura, pow (10, MULTIPLIER *  vaf_Dura ) - 1, vaf_Tumor, df.iloc[k]["gene"] )
        size = [ 100, 150, 150, 150 ]
        ax.scatter ( pow (10, MULTIPLIER *  vaf_Dura ) - 1, vaf_Tumor,  alpha = 0.7, s = size[i], 
                        color = colorlist [ df.iloc[k]["cluster_id"]], 
                        marker = shapelist[i], linewidths=0,
                        label = "{}".format ( df.iloc [k]["cluster_id"]) )
        if df.iloc[k]["gene"] in ["NF2", "AKT1", "KLF4", "TRAF7"]:
            ax.text ( pow (10, MULTIPLIER *  vaf_Dura ) - 1, vaf_Tumor, df.iloc[k]["gene"],   ha = "left", va = "bottom", fontdict = {"fontsize": 11, "fontweight" : "bold", "fontstyle": "italic"} )
        vaf_Tumor_list.append (vaf_Tumor)
        vaf_Dura_list.append (vaf_Dura)
        set_cluster_id.add ( df.iloc[k]["cluster_id"]  )

        k = k2
    
    #legend_without_duplicate_labels ( ax )

    ax.set_xticks ( [ pow (10, MULTIPLIER *  i) - 1 for i in xtick ]) 
    ax.set_xticklabels ( xtick ) 
    
    ax.set_xlim([0,  0.5 ])
    ax.set_ylim([0,  0.7])

    ax.text( ax.get_xlim()[1] / 2 , ax.get_ylim()[1] / 1.2, "Purity of Tumor = {}\nPurity of Dura = {}".format(df[df['sample_id'].str.contains('Tumor')].iloc[0]["tumour_content"], df[df['sample_id'].str.contains('Dura')].iloc[0]["tumour_content"] ), ha='center', fontsize = 8 )
    for i in range ( len(df_count[0]) ) :
        ax.text( ax.get_xlim()[1] / 2, ax.get_ylim()[1] / 1.2 - 0.03* (i + 1), "cluster{} = {}".format( df_count[0][i], int(df_count[1][i] / 2)  ), ha = 'center', fontsize = 8 )


    print ("\n\n")

    
    return df, ax










if __name__ == "__main__":
    import pandas as pd
    import matplotlib.pyplot as plt
    import palettable, argparse
    import numpy as np
    import seaborn as sns


    parser = argparse.ArgumentParser( description='The below is usage direction.')
    parser.add_argument('--Sample_ID', type=str, default="230405_2")
    parser.add_argument('--SEQUENZA_TO_PYCLONEVI_MATRIX_PATH', type=str, default="/data/project/Meningioma/31.Clonality/01.make_matrix/230405_2/230405_2.sequenza_to_pyclonevi.tsv")
    parser.add_argument('--SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH', type=str, default="/data/project/Meningioma/31.Clonality/02.pyclonevi/230405_2/230405_2.sequenza_to_pyclonevi.tsv")
    parser.add_argument('--FACETCNV_TO_PYCLONEVI_MATRIX_PATH', type=str, default="/data/project/Meningioma/31.Clonality/01.make_matrix/230405_2/230405_2.facetcnv_to_pyclonevi.tsv")
    parser.add_argument('--FACETCNV_TO_PYCLONEVI_OUTPUT_PATH', type=str, default="/data/project/Meningioma/31.Clonality/02.pyclonevi/230405_2/230405_2.facetcnv_to_pyclonevi.tsv")
    parser.add_argument('--OUTPUT_PATH_SHARED', type=str, default="/data/project/Meningioma/31.Clonality/02.pyclonevi/230405_2/230405_2.decomposed.pdf")
    parser.add_argument('--OUTPUT_DIR1', type=str, default="/data/project/Meningioma/31.Clonality/02.pyclonevi/decomposed")
    parser.add_argument('--OUTPUT_DIR2', type=str, default="/data/project/Meningioma/31.Clonality/02.pyclonevi/scaled")



    args = parser.parse_args()

    Sample_ID = args.Sample_ID
    SEQUENZA_TO_PYCLONEVI_MATRIX_PATH = args.SEQUENZA_TO_PYCLONEVI_MATRIX_PATH
    SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH = args.SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH
    FACETCNV_TO_PYCLONEVI_MATRIX_PATH = args.FACETCNV_TO_PYCLONEVI_MATRIX_PATH
    FACETCNV_TO_PYCLONEVI_OUTPUT_PATH = args.FACETCNV_TO_PYCLONEVI_OUTPUT_PATH
    OUTPUT_PATH_SHARED = args.OUTPUT_PATH_SHARED
    OUTPUT_DIR1 = args.OUTPUT_DIR1
    OUTPUT_DIR2 = args.OUTPUT_DIR2


    # df_seq_to_pycl = pd.read_csv (SEQUENZA_TO_PYCLONEVI_OUTPUT_PATH, sep = "\t")
    # df_seq_to_pycl = df_seq_to_pycl.drop_duplicates (['mutation_id'], keep = 'first').sort_values ( ['mutation_id'], axis = 0, ascending = True)[ [ "mutation_id", "cluster_id"] ].reset_index().drop ('index', axis = 1)
    #df_seq_matrix = pd.read_csv (SEQUENZA_TO_PYCLONEVI_MATRIX_PATH, sep = "\t")
     
    # # Sequenza visualization
    # fig, ax = plt.subplots( figsize=(3, 3), nrows = 1, ncols = 1 )
    # fig.subplots_adjust (wspace = 0.15, hspace = 0.1, bottom = 0.10, top = 0.92, left = 0.15, right = 0.98)
    # df_seq_integrated = pd.merge (df_seq_matrix, df_seq_to_pycl, left_on = "mutation_id", right_on = "mutation_id")
    # print ("\n\n## Decomposed - Sequenza")
    # df, ax = visualization_decomposition ( df_seq_integrated, "{} - Sequenza".format(Sample_ID), ax )
    # fig.savefig ( OUTPUT_PATH_SHARED.replace (".pdf", "-sequenza.pdf") )
    # df.to_csv ( OUTPUT_PATH_SHARED.replace (".pdf", "-sequenza.tsv"), sep = "\t", index = False)



    df_facet_to_pycl = pd.read_csv (FACETCNV_TO_PYCLONEVI_OUTPUT_PATH, sep = "\t")
    df_facet_to_pycl = df_facet_to_pycl.drop_duplicates (['mutation_id'], keep = 'first').sort_values ( ['mutation_id'], axis = 0, ascending = True)[ [ "mutation_id", "cluster_id"] ].reset_index().drop ('index', axis = 1)
    df_facet_matrix = pd.read_csv (FACETCNV_TO_PYCLONEVI_MATRIX_PATH, sep = "\t")

    # FacetCNV visualization
    fig, ax = plt.subplots( figsize=(2.5, 2.5), nrows = 1, ncols = 1 )
    fig.subplots_adjust (wspace = 0.15, hspace = 0.1, bottom = 0.12, top = 0.90, left = 0.15, right = 0.96)
    df_facet_integrated = pd.merge (df_facet_matrix, df_facet_to_pycl, left_on = "mutation_id", right_on = "mutation_id")
    print ("\n\n## Decomposed - FacetCNV")
    df, ax = visualization_decomposition ( df_facet_integrated, "{}".format(Sample_ID), ax )
    fig.savefig ( OUTPUT_PATH_SHARED.replace (".pdf", "-facetcnv.pdf"), dpi = 300)
    fig.savefig ( OUTPUT_DIR1 + "/" + str(Sample_ID) + ".decomposed-facetcnv.pdf", dpi = 300 )
    df.to_csv ( OUTPUT_PATH_SHARED.replace (".pdf", "-facetcnv.tsv"), sep = "\t", index = False)


    # FacetCNV visualization
    fig, ax = plt.subplots( figsize=(2.5, 2.5), nrows = 1, ncols = 1 )
    fig.subplots_adjust (wspace = 0.15, hspace = 0.1, bottom = 0.12, top = 0.90, left = 0.15, right = 0.96)
    df_facet_integrated = pd.merge (df_facet_matrix, df_facet_to_pycl, left_on = "mutation_id", right_on = "mutation_id")
    print ("\n\n## Decomposed & Scaled FacetCNV")
    df, ax = visualization_decomposition_scaled ( df_facet_integrated, "{}".format(Sample_ID), ax )
    fig.savefig ( OUTPUT_PATH_SHARED.replace (".pdf", "-scaled.facetcnv.pdf"), dpi = 300 )
    fig.savefig ( OUTPUT_DIR2 + "/" + str(Sample_ID) + ".scaled-facetcnv.pdf", dpi = 300 )
