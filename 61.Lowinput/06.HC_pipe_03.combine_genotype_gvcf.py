import os, subprocess, argparse

kwargs = {}

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument('--REF', type=str)
parser.add_argument('--dbSNP', type=str)
parser.add_argument('--HC_VCF_LIST', type=str)
parser.add_argument('--COMBINED_GVCF', type=str)
parser.add_argument('--GENOMICSDB_WORKSPACE', type=str)
parser.add_argument('--GENOTYPE_GVCF', type=str)
parser.add_argument('--SAMPLE_ID', type=str)
parser.add_argument('--SCRIPT_DIR', type=str)
parser.add_argument('--LOGPATH', type=str)
parser.add_argument('--TMP_DIR', type=str)
parser.add_argument('--INTERVAL', type=str)

args = parser.parse_args()

kwargs["REF"] = args.REF
kwargs["dbSNP"] = args.dbSNP
kwargs["HC_VCF_LIST"] = args.HC_VCF_LIST.split(",")
kwargs["COMBINED_GVCF"] = args.COMBINED_GVCF
kwargs["GENOMICSDB_WORKSPACE"] = args.GENOMICSDB_WORKSPACE
kwargs["GENOTYPE_GVCF"] = args.GENOTYPE_GVCF
kwargs["SAMPLE_ID"] = args.SAMPLE_ID
kwargs["SCRIPT_DIR"] = args.SCRIPT_DIR
kwargs["LOGPATH"] = args.LOGPATH
kwargs["TMP_DIR"] = args.TMP_DIR
kwargs["INTERVAL"] = args.INTERVAL

# CombineGVCFs 
# SCRIPT_combine =  "\"gatk CombineGVCFs  -R {} -O {}".format (kwargs["REF"], kwargs["COMBINED_GVCF"])
# for input_vcf in kwargs["HC_VCF_LIST"]:
#     SCRIPT_combine = SCRIPT_combine + " -V " + input_vcf
# SCRIPT_combine += "\""

# GenomicsDBImport
os.system ("rm -rf {}".format (  kwargs["GENOMICSDB_WORKSPACE"]  ) )
SCRIPT_combine =  "\"gatk  GenomicsDBImport --genomicsdb-workspace-path {} --tmp-dir {} -L {} ".format ( kwargs["GENOMICSDB_WORKSPACE"] , kwargs["TMP_DIR"], kwargs["INTERVAL"])
for input_vcf in kwargs["HC_VCF_LIST"]:
    SCRIPT_combine = SCRIPT_combine + " -V " + input_vcf
SCRIPT_combine += "\""


# GenotypeGVCFs
#SCRIPT_genotype =  "\"gatk GenotypeGVCFs  -R {} -V {} -O {}".format (kwargs["REF"], kwargs["COMBINED_GVCF"], kwargs["GENOTYPE_GVCF"] )
SCRIPT_genotype =  "\"gatk GenotypeGVCFs  -R {} -V gendb://{} -O {}".format (kwargs["REF"], kwargs["GENOMICSDB_WORKSPACE"], kwargs["GENOTYPE_GVCF"] )
SCRIPT_genotype += "\""

command_line = " ".join([ "qsub -pe smp 5 -e", kwargs["LOGPATH"], "-o", kwargs["LOGPATH"], "-N HC3_{}".format (kwargs["SAMPLE_ID"]), 
                    kwargs["SCRIPT_DIR"] + "/06.HC_pipe_03.combine_genotype_gvcf.sh",
                    "--SCRIPT_COMBINE", SCRIPT_combine, 
                    "--SCRIPT_GENOTYPE", SCRIPT_genotype,
                    "--GENOTYPE_GVCF", kwargs["GENOTYPE_GVCF"],
                    "--REF", kwargs["REF"],
                    "--dbSNP", kwargs["dbSNP"]
                        ])


os.system (command_line)


# GENOTYPE_VCF = kwargs["GENOTYPE_GVCF"].replace(".gz", "")

# grep_hash_command = "grep '^#' {} > {}".format (  GENOTYPE_VCF, GENOTYPE_VCF + ".sorted" )
# grep_sort_command = "grep -v '^#' {} | sort -k1,1V -k2n >> {}".format (GENOTYPE_VCF, GENOTYPE_VCF + ".sorted"  ) 

# # Combine the commands
# SCRIPT_sort = f"{grep_hash_command} && {grep_sort_command}"

# print (SCRIPT_sort)
# os.system ( SCRIPT_sort )
# os.system ( "mv {} {}".format ( GENOTYPE_VCF + ".sorted" ,  GENOTYPE_VCF))
# os.system ( "bgzip -c -f {} > {}".format ( GENOTYPE_VCF, kwargs["GENOTYPE_GVCF"] ) )
# os.system ( "tabix -p vcf {}".format ( kwargs["GENOTYPE_GVCF"] )  )