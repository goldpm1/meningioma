def MatrixFormation ( **kwargs ):
    import pandas as pd
    import vcf, os

    vcf_reader = vcf.Reader (open( kwargs["INPUT_VCF"], "r") )
    SigProfMatrix = pd.DataFrame ( columns = ["Project", "Sample", "ID", "Genome", "mut_type", "chrom",	"pos_start","pos_end",	"ref",	"alt"	,"Type"]) 

    for i, line in enumerate (vcf_reader):
        if len (line.REF) == len (line.ALT[0]):
            matrix = [ "Meningioma", kwargs["Sample"], line.CHROM + "_" + str(line.POS), "GRCh38", "SNP", line.CHROM.replace("chr", ""), str(line.POS), str(line.POS), line.REF, line.ALT[0], "SOMATIC "]
        elif len (line.REF) < len (line.ALT[0]):
            matrix = [ "Meningioma", kwargs["Sample"], line.CHROM + "_" + str(line.POS), "GRCh38", "INS", line.CHROM.replace("chr", ""), str(line.POS), str(line.POS), line.REF, line.ALT[0], "SOMATIC "]
        elif len (line.REF) > len (line.ALT[0]):
            matrix = [ "Meningioma", kwargs["Sample"], line.CHROM + "_" + str(line.POS), "GRCh38", "DEL", line.CHROM.replace("chr", ""), str(line.POS), str(line.POS), line.REF, line.ALT[0], "SOMATIC "]
        SigProfMatrix.loc[len(SigProfMatrix.index)] = matrix

    return SigProfMatrix



if __name__ == "__main__":
    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
    import argparse, os, glob
    import pandas as pd

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument("--RUN", type = str, default = "BY_TISSUE")
    parser.add_argument("--TISSUE", type = str, default = "")
    parser.add_argument("--SIGPROFILER_INPUT_VCF_DIR", type = str, default = "" )
    parser.add_argument("--SIGPROFILER_INPUT_MATRIX_DIR", type = str, default = "" )

    kwargs = {}
    args = parser.parse_args()

    kwargs["RUN"] = args.RUN
    kwargs["TISSUE"] = args.TISSUE
    kwargs["SIGPROFILER_INPUT_VCF_DIR"] = args.SIGPROFILER_INPUT_VCF_DIR
    kwargs["SIGPROFILER_INPUT_MATRIX_DIR"] = args.SIGPROFILER_INPUT_MATRIX_DIR

    SigProfMatrix_Total = pd.DataFrame (  columns = ["Project", "Sample", "ID", "Genome", "mut_type", "chrom",	"pos_start","pos_end",	"ref",	"alt"	,"Type"] ) 


    SIGPROFILER_INPUT_VCF_DIR_LIST = sorted (  glob.glob(kwargs["SIGPROFILER_INPUT_VCF_DIR"] + "/*.vcf") )
    for INPUT_VCF in SIGPROFILER_INPUT_VCF_DIR_LIST:
        kwargs["INPUT_VCF"] = INPUT_VCF
        if kwargs["RUN"] == "BY_TISSUE":
            kwargs["Sample"] = kwargs["TISSUE"]
        elif kwargs["RUN"] == "BY_SAMPLE":
            kwargs["Sample"] = INPUT_VCF.split(".")[0]

        SigProfMatrix = MatrixFormation ( **kwargs )
        SigProfMatrix_Total = pd.concat ( [SigProfMatrix_Total, SigProfMatrix], axis = 0 )

        if kwargs["RUN"] == "BY_SAMPLE":
            SigProfMatrix_Total.to_csv ( kwargs["SIGPROFILER_INPUT_MATRIX_DIR"] + "/" + INPUT_VCF.split("/")[-1].replace("vcf", "txt"), sep = "\t",  index = False)

    
    SigProfMatrix_Total['pos_start'] = SigProfMatrix_Total['pos_start'].astype(int)
    SigProfMatrix_Total['CHR'] = pd.Categorical(SigProfMatrix_Total['chrom'], 
                                                                        categories = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"], 
                                                                        ordered=True)
    SigProfMatrix_Total = SigProfMatrix_Total.sort_values (by = ['CHR', 'pos_start'], axis = 0).drop ("CHR", axis = 1)

    print (SigProfMatrix_Total)

    if kwargs["RUN"] == "BY_TISSUE":
        SigProfMatrix_Total.to_csv ( kwargs["SIGPROFILER_INPUT_MATRIX_DIR"] + "/" + kwargs["TISSUE"] + ".txt", sep = "\t",  index = False)