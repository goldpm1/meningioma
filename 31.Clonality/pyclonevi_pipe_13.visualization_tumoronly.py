def visualization_decomposition_1D( df, OUTPUT_SUPTITLE, ax ):
    ax.set_title(OUTPUT_SUPTITLE, fontsize = 11, fontweight='bold')
    #ax.text(0.5, 0.9, "Purity of Tumor = {}\nPurity of Dura = {}".format(df[df['sample_id'].str.contains('Tumor')].iloc[0]["tumour_content"], df[df['sample_id'].str.contains('Dura')].iloc[0]["tumour_content"] ), ha='center', fontsize = 14 )
    ax.set_xlabel("VAF_Tumor", fontdict = {"fontsize" : 8}, labelpad = 2 )
    ax.set_ylabel("Density", fontdict = {"fontsize" : 8}, labelpad = 2 )
    ax.tick_params( axis = 'x', labelsize = 8, pad = -1 )
    ax.tick_params( axis = 'y', labelsize = 8, pad = -1 )
    sns.set_style("white") 
    for axis in ['left', 'right', 'top', 'bottom']:
        ax.spines[axis].set_linewidth ( 1.25 )


    vaf_Tumor_list = []

    tabl = palettable.tableau.Tableau_20.mpl_colors
    vivid_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    if len ( set ( df["cluster_id"] ) ) >= 4:
        colorlist = ["#7288C6",  "#BCA97B", "#8D493A", "#6AEF6C" ] + [i for i in vivid_10]
    elif  len ( set ( df["cluster_id"] ) ) == 3:
        colorlist = [ "#7288C6",  "#BCA97B", "#6AEF6C"]
    elif  len ( set ( df["cluster_id"] ) ) == 2:
        colorlist = [ "#7288C6", "#BCA97B",]
    else:
        colorlist = ["#7288C6]
    shapelist = ["o", "s", "^", "v", ">", "<"]

    
    df[ ["CHR", "POS", "ALT", "REF"] ] = df["mutation_id"].str.split("_", expand = True)
    df["CHR"] = pd.Categorical(df["CHR"], 
                                                categories = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"], 
                                                ordered=True)
    df = df.sort_values (by = ['cluster_id', 'CHR', 'POS'], axis = 0).drop (["CHR", "POS"], axis = 1).reset_index(drop = True)
    df["VAF"] = round(df["alt_counts"].astype(float) / (df["ref_counts"].astype(float) + df["alt_counts"].astype(float)), 3)

    
    # df["cluster_id"]를 세어서 count 해주는 dataframe으로 만들어줍니다.
    df_count = df["cluster_id"].value_counts().sort_index().reset_index()
    df_count.columns = ["cluster_id", "count"]
    print(df_count)


    
    # cluster_id별로 VAF에 대한 density plot을 그립니다.
    for idx, cluster in enumerate ( sorted(df["cluster_id"].unique()) ):
        cluster_vaf = df[df["cluster_id"] == cluster]["VAF"]
        sns.kdeplot ( cluster_vaf,   ax = ax,  label=f"cluster {cluster}",  color=colorlist[idx % len(colorlist)],  fill=True,  alpha=0.3, linewidth=1.5  )
        ax.scatter(cluster_vaf, [0]*len(cluster_vaf), color=colorlist[idx % len(colorlist)], marker=shapelist[idx % len(shapelist)], s = 100, alpha=0.7, edgecolors='black', linewidths=0.5)

    
    # "gene" 컬럼이 ["NF2", "AKT1", "KLF4", "TRAF7"]에 해당하는 row의 index를 찾아서, 해당 VAF 위치에 gene 이름을 표시
    for k in df.index:
        if df.iloc[k]["gene"] in ["NF2", "AKT1", "KLF4", "TRAF7"]:
            x = df.iloc[k]["VAF"]
            y = ax.get_ylim()[1] / 5
            ax.text(x, y, df.iloc[k]["gene"], ha="center", va="bottom", fontdict={"fontsize": 11, "fontweight": "bold", "fontstyle": "italic"})
            # y1 = ax.get_ylim(1)/6 부터 y2=0까지 화살표를 그려줘
            ax.annotate( '',  xy=(x, 0),   xytext=(x, ax.get_ylim()[1] / 5),   arrowprops=dict(arrowstyle='->', color='black', lw=1.5)        )






    if ("190426" in OUTPUT_SUPTITLE) | ("221102" in OUTPUT_SUPTITLE) | ("230127" in OUTPUT_SUPTITLE) | ("230405" in OUTPUT_SUPTITLE) | ("230419" in OUTPUT_SUPTITLE) | ("230802" in OUTPUT_SUPTITLE) | ("230822" in OUTPUT_SUPTITLE) | ("230920" in OUTPUT_SUPTITLE):
        ax.set_xlim([0,  0.7]);  
    elif "220930" in OUTPUT_SUPTITLE:
        ax.set_xlim([0,  0.55]); 
    else:
        ax.set_xlim([0,  0.45]); 

    ax.text( ax.get_xlim()[1] / 2 , ax.get_ylim()[1] / 1.2, "Purity of Tumor = {}".format(df[df['sample_id'].str.contains('Tumor')].iloc[0]["tumour_content"] ), ha='center', fontsize = 8 )
    for i in range ( df_count.shape[0] ) :
        ax.text( ax.get_xlim()[1] / 2, ax.get_ylim()[1] / 1.2 - (ax.get_ylim()[1]/10)* (i + 1), "cluster{} = {}".format( df_count.iloc[i]["cluster_id"],  df_count.iloc[i]["count"]  ), ha = 'center', fontsize = 8 )


    print ("\n\n")
    
    return df, ax
        








if __name__ == "__main__":
    import pandas as pd
    import matplotlib.pyplot as plt
    import palettable, argparse
    import numpy as np
    import seaborn as sns
    import matplotlib as mpl
    import os, glob

    # List all available fonts
    file_paths = glob.glob( "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/arial*" )
    absolute_file_paths = [os.path.abspath(file_path) for file_path in file_paths]
    for absolute_file_path in absolute_file_paths:
        mpl.font_manager.fontManager.addfont( absolute_file_path )
    plt.rcParams["font.family"] = 'Arial'


    parser = argparse.ArgumentParser( description='The below is usage direction.')
    parser.add_argument('--Sample_ID', type=str, default="230405_2")
    parser.add_argument('--FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH', type=str, default="/data/project/Meningioma/31.Clonality/01.make_matrix/230405_2/230405_2.facetcnv_to_pyclonevi.tsv")
    parser.add_argument('--FACETCNV_TO_PYCLONEVI_OUTPUT_PATH', type=str, default="/data/project/Meningioma/31.Clonality/02.pyclonevi/230405_2/230405_2.facetcnv_to_pyclonevi.tsv")
    parser.add_argument('--OUTPUT_VIS_FIG1', type=str, default="")
    parser.add_argument('--OUTPUT_VIS_FIG2', type=str, default="")
    parser.add_argument('--OUTPUT_VIS_DF', type=str, default="")
    parser.add_argument('--OUTPUT_DIR1', type=str, default="/data/project/Meningioma/31.Clonality/02.pyclonevi/decomposed")
    parser.add_argument('--OUTPUT_DIR2', type=str, default="/data/project/Meningioma/31.Clonality/02.pyclonevi/scaled")



    args = parser.parse_args()

    Sample_ID = args.Sample_ID
    FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH = args.FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH
    FACETCNV_TO_PYCLONEVI_OUTPUT_PATH = args.FACETCNV_TO_PYCLONEVI_OUTPUT_PATH
    OUTPUT_VIS_FIG1 = args.OUTPUT_VIS_FIG1
    OUTPUT_VIS_FIG2 = args.OUTPUT_VIS_FIG2
    OUTPUT_VIS_DF = args.OUTPUT_VIS_DF
    OUTPUT_DIR1 = args.OUTPUT_DIR1
    OUTPUT_DIR2 = args.OUTPUT_DIR2




    df_facet_to_pycl = pd.read_csv (FACETCNV_TO_PYCLONEVI_OUTPUT_PATH, sep = "\t")
    df_facet_to_pycl = df_facet_to_pycl.drop_duplicates (['mutation_id'], keep = 'first') [ [ "mutation_id", "cluster_id"] ].reset_index().drop ('index', axis = 1)
    df_facet_matrix = pd.read_csv (FACETCNV_TO_PYCLONEVI_1D_MATRIX_PATH, sep = "\t")   

    # FacetCNV visualization
    fig, ax = plt.subplots( figsize=(2, 2), nrows = 1, ncols = 1 )
    fig.subplots_adjust (wspace = 0.15, hspace = 0.1, bottom = 0.14, top = 0.88, left = 0.18, right = 0.96)
    df_facet_integrated = pd.merge (df_facet_matrix, df_facet_to_pycl, left_on = "mutation_id", right_on = "mutation_id")  # merge 해서 cluster_id 추가
    print ("\n\n## Decomposed - FacetCNV")
    df, ax = visualization_decomposition_1D ( df_facet_integrated, "{}".format(Sample_ID), ax )
    fig.savefig ( OUTPUT_VIS_FIG1, dpi = 300)
    fig.savefig ( OUTPUT_DIR1 + "/" + str(Sample_ID) + ".facetcnv.pdf", dpi = 300 )
    print (OUTPUT_VIS_DF)
    df.to_csv ( OUTPUT_VIS_DF, sep = "\t", index = False)

