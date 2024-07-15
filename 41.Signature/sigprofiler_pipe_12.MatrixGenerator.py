from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import argparse

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument("--PROJECT", type = str, default = "Meningioma")
parser.add_argument("--OUTPUT_DIR", type = str, default = "/data/project/Meningioma/41.Signature/01.SigProfiler/11.matrix/Shared")

kwargs = {}
args = parser.parse_args()

kwargs["PROJECT"] = args.PROJECT
kwargs["OUTPUT_DIR"] = args.OUTPUT_DIR

matrices = matGen.SigProfilerMatrixGeneratorFunc ( kwargs["PROJECT"],  # project
                                                                                    "GRCh38",   # reference_genome
                                                                                    kwargs["OUTPUT_DIR"],  # path_to_input_files 
                                                                                    exome=False, 
                                                                                    bed_file=None, 
                                                                                    chrom_based=False, 
                                                                                    plot=True, 
                                                                                    tsb_stat=False, 
                                                                                    seqInfo=True)
                                                                                    