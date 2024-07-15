import pandas as pd
import numpy as np
from os import path
import argparse
from collections import Counter

pd.set_option('display.max_seq_items', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 500)

parser = argparse.ArgumentParser(description='Here is usage direction.')
parser.add_argument(    '--INPUTVCF', default="/home/goldpm1/Meningioma/04.mutect/04.rescue/221102_Tumor.MT2.FMC.HF.RMBLACK.vep.rescue.vcf")
parser.add_argument(    '--OUTPUT_MAF', default="/home/goldpm1/Meningioma/04.mutect/06.maf/01.shared_yes/221102_Tumor.MT2.FMC.HF.RMBLACK.vep.vcf")
parser.add_argument(    '--SELECTED_DB', default="")

args = parser.parse_args()
INPUTVCF = args.INPUTVCF
OUTPUT_MAF = args.OUTPUT_MAF
SELECTED_DB = args.SELECTED_DB
SELECTED_DB = SELECTED_DB.split(",")


input_file = open(INPUTVCF, "r")


ID_list = []
ID_DA = []       # Depth, Alt를 tuple in list 형태로 저장
ID_DP = []       # Variant가 있는 경우 → Depth의 list in list
ID_Alt = []      # Variant가 있는 경우 → Alt의 list in list
ID_VAF = []      # Variant가 있는 경우 → VAF의 list in list

VEP_format = []

matrix = []
maf_colnames = ['Chromosome', 'Start_Position', 'End_Position',
                'Tumor_Sample_Barcode', "ID", 'Hugo_Symbol', 'Reference_Allele',
                'Tumor_Seq_Allele2',  "Alt", "Depth", "BIOTYPE", 'Variant_Classification', "Impact",
                'tx', 'exon', "rsID", 'txChange', 'aaChange', 'Variant_Type', 'sample_id',
                'hgnc_symbol', 'Entrez', 'ens_id', 'Entrez_Gene_Id',
                "clinvar_MedGen_id", "clinvar_OMIM_id", "clinvar_Orphanet_id", "clinvar_clnsig", "clinvar_hgvs", "clinvar_id", "clinvar_review", "clinvar_trait", "clinvar_var_source",
                "CLIN_SIG", "SIFT", "PolyPhen",
                "SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred",
                "ClinPred_pred", "ClinPred_rankscore", "ClinPred_score", "DANN_rankscore", "DANN_score",
                "CADD", "GERP++_RS", "FATHMM_pred", "FATHMM_score", 'MetaSVM_pred',  "MetaSVM_rankscore", "MetaSVM_score", "MutationTaster_pred", 'MutationAssessor_pred',
                "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "SpliceAI_score",
                'MOTIF_NAME', 'MOTIF_POS', "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS", 'FunMotifs', "Nearest_gene",
                "KRG", "K1", "gnomAD", "dbSNP"]      # vep 돌리면 나오는 68개 
# noncoding_colnames = ["ChIP", "DNase", "PWM", "Footprint", "QTL", "PWM_matched", "Footprint_matched", "ranking_probability", "Non_Coding_Score", "Non_Coding_Groups", "Coding_Score", "Coding_Group", "EA_enhancer", "GH_promoter", "GH_enhancer", "RefSeq_promoter", "GH_promoter_enhancer", "ENSEMBL_promoter", "mTL_miRNA",
#                       "greendb_id", "greendb_stdtype", "greendb_dbsource", "greendb_genes", "green_constraint", "greendb_level", "ANN"] + SELECTED_DB  # 26개 + alpha
noncoding_colnames = []
maf_colnames = maf_colnames + noncoding_colnames

maf_Variant_Classification = {"missense_variant": "Missense_Mutation", "stop_gained": "Nonsense_Mutation", "synonymous_variant": "Synonymous_Mutation"}
maf_Variant_Type = {"insertion": "INS", "deletion": "DEL", "SNV": "SNP"}
df_maf = pd.DataFrame(columns=maf_colnames)


def GenotypeParsing(line, i, ID_list):
    global ID_DA
    global ID_DP
    global ID_Alt
    global ID_VAF
    sample_variantNo = []
    ID_VAF_line = []

    formatdata = line[8].split(':')
    sampledata = {}
    for j in range(len(formatdata)):
        sampledata[formatdata[j]] = line[i].split(':')[j]

    if sampledata['GT'] in ["0", "1", "0/1", "0|1", "1/0", "1|0", "0/2", "0|2", "1/1", "1|1", "2/2", "2|2"]:
        try:
            ref = int(sampledata['AD'].split(',')[0])
        except:
            ref = 0

        try:
            alt = int(sampledata['AD'].split(',')[1])
        except:
            alt = 0
        dp = ref + alt

        sample_variantNo =  ID_list[i - 9]

        if alt == 0:
            return sample_variantNo, alt, dp

        ID_DA[i-9].append((alt, dp))
        ID_DP[i-9].append(dp)
        ID_Alt[i-9].append(alt)
        ID_VAF[i-9].append(int(alt) / dp)
        ID_VAF_line.append(round(float(int(alt) / dp), 2))

        return sample_variantNo, alt, dp
    else:   # 매우 이상한 경우
        return [0], 0, 0


def Decision_meaningful_noncoding(info_dict, info_line):

    if "ANN" in info_dict:      # GreenDB annotation이 붙었을 경우.  GreenDB는 절대적으로 신뢰
        if "," in info_dict["greendb_genes"]:
            a = info_dict["greendb_genes"].split(",")
        else:
            a = [info_dict["greendb_genes"]]

        for aa in info_dict["greendb_stdtype"].split(","):
            if ("promoter" in aa):
                return True, "promoter", a
                # [ANN_line.split("|")[3]]          # GreenDB의 genename을 같이 반환
        for aa in info_dict["greendb_stdtype"].split(","):
            if ("enhancer" in aa):
                return True, "enhancer", a          # GreenDB의 genename을 같이 반환

    check_db = {"RegulomeDB": 0, "Fathmm_MKL": 0,
                "GH+Ensembl": 0, "regBase": 0}

    # GH + Ensembl
    g1 = ""
    consequence = ""
    for GH_Ensembl in ["EA_enhancer", "GH_promoter", "GH_enhancer", "RefSeq_promoter", "GH_promoter_enhancer", "ENSEMBL_promoter"]:
        if info_dict[GH_Ensembl] != ".":       # 뭐라도 있으면 1점 추가.  그리도 g1 list에 추가해서 반환
            check_db["GH+Ensembl"] = check_db["GH+Ensembl"] + 1
            if "promoter" in GH_Ensembl:
                consequence = "promoter"
            if (consequence == "") & ("enhancer" in GH_Ensembl):
                consequence = "enhancer"

            # genelist는 다 합쳐주고 나중에 set로 만들 예정
            g1 = g1 + "," + info_dict[GH_Ensembl]

    if consequence != "":            # 뭐 하나라도 걸리면 돌려보내준다
        g1 = list(set(g1.split(",")))
        g1 = list(filter(None, g1))
        return True, consequence, g1

    # RegulomeDB
    for keys in ["ChIP", "DNase", "PWM", "Footprint", "QTL", "PWM_matched", "Footprint_matched"]:
        if keys in info_dict:
            if info_dict[keys] == "True":
                check_db["RegulomeDB"] = check_db["RegulomeDB"] + 0.2

    if info_dict["ranking_probability"] != ".":
        # 의미 있다고 생각하고 돌려보내긴 하는데 gene 이름을 주기가 힘들다.  nearest gene이 있으면 좋겠다
        return True, "noncoding", [info_dict["Nearest_gene"]]

    # Fathmm_MKL
    if (info_dict["Non_Coding_Score"] != ".") & (info_dict["Non_Coding_Score"] != "Non_Coding_Score"):
        try:
            t = float(info_dict["Non_Coding_Score"])
        except:
            print(info_dict, info_line)
        # 의미 있다고 생각하고 돌려보내긴 하는데 gene 이름을 주기가 힘들다.  nearest gene이 있으면 좋겠다
        if float(info_dict["Non_Coding_Score"]) > 0.5:
            return True, "noncoding", [info_dict["Nearest_gene"]]

    # regBase (DeepSEA, ReMM, PAFA)
    for DB in ["DeepSEA", "ReMM", "PAFA"]:
        if info_dict[DB + "_PHRED"] != ".":
            # 일단 이정도로 돌려보내긴 하는데 gene이름을 찾기가 어렵다. nearest gene이 있으면 좋겠다
            if float(info_dict[DB + "_PHRED"]) > 5:
                return True, "noncoding", [info_dict["Nearest_gene"]]

    return False, "", ""


def Multiple_Consequence(li):    # 여러 개 annotation이 붙은 경우 어떤 것을 우선순위로 줄지
    for i in li:
        if ("splice_donor" in i) | ("splice_acceptor" in i):
            return "Splice_site"
        elif "frameshift_variant" in i:
            return "frameshift_variant"
        elif "stop_gained" in i:
            return "stop_gained"
        elif "stop_lost" in i:
            return "stop_lost"
        elif "missense" in i:
            return "missense_variant"

    return "multiple_meaningless"


Variant_Classification_set = set([])

cnt = 0
while True:
    line = input_file.readline()
    # line = line.decode('utf-8')             # input이 gzip일 경우
    line = line.rstrip('\n')

    if len(line) < 2:
        break

    if line[0:14] == "##INFO=<ID=CSQ":
        VEP_format = line.split("Format: ")[1].split("\"")[0].split("|")

    if line[0:2] == "##":                      # Header만 따로 저장이 필요할 경우
        continue

    original_line = line
    line = line.split()

    if line[0] == "#CHROM":
        for i in range(9, len(line)):
            ID_list.append(line[i])
            ID_DA.append([])
            ID_DP.append([])
            ID_Alt.append([])
            ID_VAF.append([])
        continue

    try:
        if 'AD' not in line[8]:
            continue
    except:
        print (line)
        break

    cnt = cnt + 1

    # 기본 VCF parsing하기
    CHR = line[0]
    POS = int(line[1])
    REF = line[3]
    ALT = line[4]


    # VEP 정보 parsing하기  (info_dict)
    info_list = line[7].split(';')
    info_dict = {}
    info_dict["Nearest_gene"] = "-"
    for key in ["KRG", "K1", "gnomAD", "dbSNP"] + noncoding_colnames:
        info_dict[key] = "."

    for i in range(len(info_list)):
        if "CSQ" in info_list[i]:
            continue

        if '=' not in info_list[i]:  # ?=?  형식으로 안 돼있는 것도 있다
            info_dict[info_list[i]] = info_list[i]
        else:
            info_dict[info_list[i].split('=')[0]] = info_list[i].split('=')[1]

    
    
    
    ######## CSQ= 다음이 vep가 annotation해준 것.  그런데 multiple traanscript일 경우 , 로 구분되어 있다 ########
    VEP_verbose = line[7].split("CSQ=")[1].split(";")[0].split(",")

    for j in range(len(VEP_verbose)):  # multiple transcript일 경우 하나씩 도는 것
        VEP_dict = {}
        for k in maf_colnames:         # 일단 다 빈칸으로 채워넣고
            VEP_dict[k] = ''
        # |로 구분된 CSQ 정보를 split 해서 채워넣는다
        
        for k in range( len(VEP_verbose[j].split('|')) ):
            #print ( "VEP_format[k] = {}\tVER_verbose[j].split("")[k] = {}".format ( VEP_format[k], VEP_verbose[j].split("|")[k] ))
            try:
                VEP_dict[ VEP_format[k] ] = VEP_verbose[j].split('|')[k]
                #print ( "VEP_format[k] = {}\tVER_verbose[j].split("")[k] = {}".format ( VEP_format[k], VEP_verbose[j].split('|')[k]))
            except:
                #print(" No | in VEP_verbose[j], k ({})".format ( k, VEP_verbose[j].split('|')[k] ) )
                continue

        # Population 정보 update
        if "dbSNP_K1" in VEP_dict.keys():
            if VEP_dict["dbSNP_K1"] != "":
                info_dict["K1"] = VEP_dict["dbSNP_K1"]
        if "dbSNP_KRG" in VEP_dict.keys():
            if VEP_dict["dbSNP_KRG"] != "":
                info_dict["KRG"] = VEP_dict["dbSNP_KRG"]
        if "dbSNP_EAS_AF" in VEP_dict.keys():
            if VEP_dict["dbSNP_EAS_AF"] != "":      # dbSNP frequency 채워주기
                info_dict["dbSNP"] = VEP_dict["dbSNP_EAS_AF"]
        if "gnomAD_AF" in VEP_dict.keys():
            if VEP_dict["gnomAD_AF"] != "":
                info_dict["gnomAD"] = VEP_dict["gnomAD_AF"]


        if (VEP_dict["BIOTYPE"] in ["protein_coding", "promoter", "CTCF_binding_site", "enhancer"]) | (VEP_dict["Consequence"] == "TF_binding_site_variant"):
            # 여기는 regulatory_region_variant -CTCF_binding_site로 되어 있다
            if VEP_dict["BIOTYPE"] == "CTCF_binding_site":
                VEP_dict["Consequence"] = "CTCF_binding_site"
            #     VEP_dict["SYMBOL"] = info_dict["Nearest_gene"]
            # if VEP_dict["Consequence"] == "TF_binding_site_variant":
            #     VEP_dict["SYMBOL"] = info_dict["Nearest_gene"]

            #print(VEP_dict["BIOTYPE"], VEP_dict["Consequence"])

            # if VEP_dict["Consequence"] in ["intron_variant", "upstream_gene_variant", "synonymous_variant", "downstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant"]:          # 의미 없는건 버리자
            #     continue

            # mutliple 인 경우에는 우선순위 높은 거로 하나만 살리자
            if "&" in VEP_dict["Consequence"]:
                VEP_dict["Consequence"] = Multiple_Consequence(VEP_dict["Consequence"].split("&"))            # 우선순위 높은 걸로 하나만 살리자

            for i in range(10, 11):        # Multisample일 수 있으니까
                sample_variantNo, alt, depth = GenotypeParsing(line, i, ID_list)
                
                if alt != 0:
                    filtered_list = [value for value in [VEP_dict["SpliceAI_pred_DS_AG"], VEP_dict["SpliceAI_pred_DS_AL"],  VEP_dict["SpliceAI_pred_DS_DG"], VEP_dict["SpliceAI_pred_DS_DL"]] if value != '' and value != None]
                    SpliceAI_score = ''
                    if len(filtered_list) != 0:
                        filtered_list = [float(i) for i in filtered_list]
                        SpliceAI_score = np.max(np.array(filtered_list))
                        if SpliceAI_score > 0.8:
                            VEP_dict["Consequence"] = "Splice_Site"

                    # VEP에서 준 VEP_dict["VARIANT_CLASS"]를 maf 양식에 맞춰 바꾸는 과정     {"insertion": "INS", "deletion": "DEL", "SNV": "SNP"}
                    if VEP_dict["VARIANT_CLASS"] in maf_Variant_Type.keys():
                        VEP_dict["VARIANT_CLASS"] = maf_Variant_Type[VEP_dict["VARIANT_CLASS"]]

                    # VEP에서 준 VEP_dict["Consequence"]를 maf 양식에 맞춰 바꾸는 과정    {"missense_variant": "Missense_Mutation","stop_gained": "Nonsense_Mutation", "regulatory_region_variant": "Translation_Start_Site", "frameshift_variant":"Frame_Shift_Ins"}
                    if VEP_dict["Consequence"] in maf_Variant_Classification.keys():
                        VEP_dict["Consequence"] = maf_Variant_Classification[VEP_dict["Consequence"]]
                    if VEP_dict["Consequence"] == "frameshift_variant":
                        if VEP_dict["VARIANT_CLASS"] == "INS":
                            VEP_dict["Consequence"] = "Frame_Shift_Ins"
                        if VEP_dict["VARIANT_CLASS"] == "DEL":
                            VEP_dict["Consequence"] = "Frame_Shift_Del"

                    Variant_Classification_set.add( VEP_dict["BIOTYPE"] + "_" + VEP_dict["Consequence"])


                    matrix_add = [CHR, POS - 1, POS, ID_list[i-9], VEP_dict["SYMBOL"] + "_" + str(POS), VEP_dict["SYMBOL"], REF, ALT, alt, depth, VEP_dict["BIOTYPE"], VEP_dict["Consequence"], VEP_dict["IMPACT"], VEP_dict["Feature"], VEP_dict["EXON"], VEP_dict["Existing_variation"].split("&")[0], VEP_dict["HGVSc"], VEP_dict["HGVSp"], VEP_dict["VARIANT_CLASS"], VEP_dict["SYMBOL"],
                                            None, None, VEP_dict["SYMBOL"], None, VEP_dict["clinvar_MedGen_id"], VEP_dict["clinvar_OMIM_id"], VEP_dict["clinvar_Orphanet_id"], VEP_dict["clinvar_clnsig"], VEP_dict["clinvar_hgvs"], VEP_dict["clinvar_id"], VEP_dict["clinvar_review"], VEP_dict["clinvar_trait"], VEP_dict["clinvar_var_source"],
                                            VEP_dict["CLIN_SIG"], VEP_dict["SIFT"], VEP_dict["PolyPhen"], VEP_dict["SIFT_pred"], VEP_dict["Polyphen2_HDIV_pred"], VEP_dict["Polyphen2_HVAR_pred"],
                                            VEP_dict["CADD_phred"], VEP_dict["GERP++_RS"], VEP_dict["FATHMM_pred"], VEP_dict["FATHMM_score"], VEP_dict["MetaSVM_pred"], VEP_dict["MetaSVM_rankscore"], VEP_dict["MetaSVM_score"], VEP_dict["MutationTaster_pred"], VEP_dict['MutationAssessor_pred'], VEP_dict["ClinPred_pred"], VEP_dict["ClinPred_rankscore"], VEP_dict["ClinPred_score"], VEP_dict["DANN_rankscore"], VEP_dict["DANN_score"],
                                            VEP_dict['SpliceAI_pred_DS_AG'], VEP_dict['SpliceAI_pred_DS_AL'], VEP_dict['SpliceAI_pred_DS_DG'], VEP_dict['SpliceAI_pred_DS_DL'],
                                            SpliceAI_score,
                                            VEP_dict['MOTIF_NAME'], VEP_dict['MOTIF_POS'], VEP_dict["MOTIF_SCORE_CHANGE"], VEP_dict["TRANSCRIPTION_FACTORS"], VEP_dict["FunMotifs"], info_dict["Nearest_gene"],
                                            info_dict["KRG"], info_dict["K1"], info_dict["gnomAD"], info_dict["dbSNP"]] + [info_dict[a] for a in noncoding_colnames]

                    matrix.append(matrix_add)



    #################  Noncoding 기준에 맞으면 Consequence == "noncoding"으로 삽입해줌 #################
    # tf, consequence, genename = Decision_meaningful_noncoding(info_dict, line[7])

    # if tf == True:
    #     sample_variantNo, alt, depth = GenotypeParsing(line, 9, ID_list)

    #     if ID_list == [0]:  # 엄청 이상한 경우
    #         continue

    #     for gene in genename:
    #         temp = ["." for i in range(len(maf_colnames))]

    #         try:
    #             temp[0:4] = [CHR, POS - 1, POS, sample_variantNo[0]]
    #         except:
    #             print("Error occurred")
    #             print(line, CHR, POS, sample_variantNo)
    #             break
    #         temp[4:6] = [gene + "_" + str(POS), gene]
    #         temp[19] = gene
    #         temp[22] = gene

    #         if len(REF) == len(ALT):
    #             temp[18] = "SNP"
    #         else:
    #             temp[18] = "Indel"
    #         temp[len(temp) - len(noncoding_colnames): len(temp)
    #              ] = [info_dict[a] for a in noncoding_colnames]

    #         temp[6: 13] = [REF, ALT, alt, depth, "noncoding", consequence,
    #                        "MODIFIER"]       # e.g)    promoter, enhancer"'
    #         matrix.append(temp)


df_maf = df_maf.append(pd.DataFrame.from_records(matrix, columns=maf_colnames))

df_maf.to_csv(OUTPUT_MAF, index=False, na_rep='.', sep='\t', header=True)

input_file.close()
