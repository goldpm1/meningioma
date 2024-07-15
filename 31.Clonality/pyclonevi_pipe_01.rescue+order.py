def order_by_chrpos ( df  ):
    # 읽고 chr_pos 로 정렬해줌
    df = df.drop_duplicates(["mutation_id", "sample_id"],  keep = 'first')
    df [["chr", "pos"]]  = df["mutation_id"].str.split("_", 1 , expand = True)  # chr, pos로 펼쳐준 다음 예쁘게 sorting 하기
    df.loc[:,"chr"] = df.loc[:,"chr"].str.replace ("chr", "")
    df = df.astype({'chr': 'int'} )

    df.sort_values ( ['chr', "pos", 'sample_id'], axis = 0, ascending = True, inplace = True)
    df = df.drop ( ["chr", "pos"], axis = 1) # 다시 그 column 없애주기
    df = df.reset_index(drop = True) 

    return df


def rescue_unique_mutation (df):
    import numpy as np
    import copy

    # 빈 데이터프레임 생성
    df_new = pd.DataFrame(columns=df.columns)

    # TUMOR, DURA의 purity를 dictionary에 담기
    sample_purity_dict = {}
    for i in pd.DataFrame ( df.groupby(["tumour_content"])["sample_id"].value_counts() ).index:
        sample_purity_dict [ i[1] ] = i[0]

    # 홀로인 mutation_id의 list
    df_cnt = df.groupby("mutation_id").count() 
    list_cnt = np.array ( df_cnt.iloc[:, 0].tolist() )
    incomplete_mutation_id_list =  df_cnt.index [np.where (list_cnt != len (sample_purity_dict.keys())) [0]  ].tolist()     # 1인 애들을 뽑아오기


    ####### 완벽하지 않은 애들은 복붙해주기 #######3
    for k in range ( df.shape[0] ):
        if ( df.iloc[k]["mutation_id"] in incomplete_mutation_id_list ):
            row_series = df.iloc[k]

            for s in sample_purity_dict.keys():    # "220930_Dura", "220930_Tumor"
                if s not in row_series["sample_id"]:  # 없다면 
                    u = row_series["sample_id"]
                    row_series_copy = copy.deepcopy (row_series)
                    row_series_copy ["sample_id"] = row_series_copy["sample_id"].replace ( u, s)   # 220930_Dura -> 220930_Tumor
                    row_series_copy ["alt_counts"] = 0
                    row_series_copy ["tumour_content"] = sample_purity_dict [ row_series_copy["sample_id"] ]
                    
                    df_new = df_new.append ( row_series_copy, ignore_index=True)

    # 합쳐주기
    df_total = df.append (df_new, ignore_index = True).sort_values ( ['mutation_id', 'sample_id'], axis = 0, ascending = True).reset_index(drop = True) 

    return df_total


if __name__ == "__main__":
    import pandas as pd
    import argparse
    pd.options.mode.chained_assignment = None  # Suppress the SettingWithCopyWarning
    import warnings
    warnings.simplefilter (action = 'ignore')

    parser = argparse.ArgumentParser( description='The below is usage direction.')
    parser.add_argument('--Sample_ID', type=str, default="220930")
    parser.add_argument('--TISSUE', type=str, default="Tumor")
    parser.add_argument('--SEQUENZA_TO_PYCLONEVI_MATRIX_PATH', type=str, default="")
    parser.add_argument('--FACETCNV_TO_PYCLONEVI_MATRIX_PATH', type=str, default="")
    parser.add_argument('--RESCUE_UNIQUEMUTATION', type=bool)

    args = parser.parse_args()

    Sample_ID = args.Sample_ID
    TISSUE = args.TISSUE
    SEQUENZA_TO_PYCLONEVI_MATRIX_PATH = args.SEQUENZA_TO_PYCLONEVI_MATRIX_PATH
    FACETCNV_TO_PYCLONEVI_MATRIX_PATH = args.FACETCNV_TO_PYCLONEVI_MATRIX_PATH
    RESCUE_UNIQUEMUTATION = bool (args.RESCUE_UNIQUEMUTATION)


    df_SEQUENZA_TO_PYCLONEVI = pd.read_csv (SEQUENZA_TO_PYCLONEVI_MATRIX_PATH, sep = "\t", names =["mutation_id", "sample_id", 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn', 'minor_cn', 'tumour_content', 'types' , 'gene', 'variant_classification'] )
    df_FACETCNV_TO_PYCLONEVI = pd.read_csv (FACETCNV_TO_PYCLONEVI_MATRIX_PATH, sep = "\t", names =["mutation_id", "sample_id", 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn', 'minor_cn', 'tumour_content', 'types' , 'gene', 'variant_classification'] )

    #살려줄거면 살려주기
    if RESCUE_UNIQUEMUTATION == True:
        df_SEQUENZA_TO_PYCLONEVI = rescue_unique_mutation ( df_SEQUENZA_TO_PYCLONEVI )
        df_FACETCNV_TO_PYCLONEVI = rescue_unique_mutation ( df_FACETCNV_TO_PYCLONEVI )

    # 예쁘게 정렬
    df_SEQUENZA_TO_PYCLONEVI = order_by_chrpos ( df_SEQUENZA_TO_PYCLONEVI )
    df_FACETCNV_TO_PYCLONEVI = order_by_chrpos ( df_FACETCNV_TO_PYCLONEVI )

    # SAVE
    df_SEQUENZA_TO_PYCLONEVI.to_csv (SEQUENZA_TO_PYCLONEVI_MATRIX_PATH, sep = "\t", index = False)
    df_FACETCNV_TO_PYCLONEVI.to_csv (FACETCNV_TO_PYCLONEVI_MATRIX_PATH, sep = "\t", index = False)