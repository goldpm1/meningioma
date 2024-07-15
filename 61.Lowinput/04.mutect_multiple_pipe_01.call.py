import os, subprocess, argparse

kwargs = {}

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument('--REF', type=str)
parser.add_argument('--BAM_DIR_LIST', type=str)
parser.add_argument('--normal', type=str)
parser.add_argument('--panel_of_normals', type=str)
parser.add_argument('--germline_resource', type=str)
parser.add_argument('--O', type=str)
parser.add_argument('--temp_dir', type=str)

args = parser.parse_args()

kwargs["REF"] = args.REF
kwargs["BAM_DIR_LIST"] = args.BAM_DIR_LIST.split(",")
kwargs["normal"] = args.normal
kwargs["panel_of_normals"] = args.panel_of_normals
kwargs["germline_resource"] = args.germline_resource
kwargs["O"] = args.O
kwargs["temp_dir"] = args.temp_dir


# Mutect Multiple
os.system ("rm -rf {}".format (  kwargs["O"]  ) )
SCRIPT =  "gatk  Mutect2 -R {} -normal {} --panel-of-normals {} --germline-resource {} -O {} --tmp-dir {} ".format ( 
    kwargs["REF"], kwargs["normal"], kwargs["panel_of_normals"], kwargs["germline_resource"], kwargs["O"], kwargs["temp_dir"]    )
for input_bam in kwargs["BAM_DIR_LIST"]:
    if input_bam != "":
        SCRIPT = SCRIPT + " -I " + input_bam
#SCRIPT += "\""

print (SCRIPT)
os.system (SCRIPT)