#!/usr/bin/env python3
import os, subprocess
import shutil
import time

CURRENT_PATH = os.getcwd()
logPath = os.path.join(CURRENT_PATH, "log")

PROJECT_DIR = "/data/project/Meningioma"
BAM_DIR = os.path.join(PROJECT_DIR, "02.Align")
MUTECT_DIR = os.path.join(PROJECT_DIR, "04.mutect")
PON = "/data/public/GATK/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
#REF = "/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa"    # contig: 196
REF = "/home/goldpm1/reference/genome.fa"    # contig: 3367
hg = "hg38"
gnomad = "/data/public/GATK/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
TMP_PATH = os.path.join(BAM_DIR, "temp")
INTERVAL = "/home/goldpm1/resources/Agilent_SureSelectXT_Human_All_Exon_Kit_V5_hg38/S04380110_Covered.bed"

if not os.path.exists(logPath):
    os.makedirs(logPath)

for sublog in ["11.multipleMT", "02.FMC_HF_RMBLACK"]:
    sublog_path = logPath + "/"  + sublog
    if os.path.exists(sublog_path):       
        shutil.rmtree(sublog_path)
    if not os.path.exists(sublog_path):
        os.makedirs(sublog_path, exist_ok=True)


######################################################################################################
def get_slurm_jobid(slurm_cmd):
    # slurm_command 문자열을 바로 실행하고, jobid를 1회만 추출
    try:
        result = subprocess.run(
            slurm_command,
            shell=True, capture_output=True, text=True,   check=True
        )
        output = result.stdout.strip()
        # 잡아이디 추출
        import re
        m = re.search(r"Submitted batch job (\d+)", output)
        if m:
            job_id = m.group(1)
        else:
            print("SLURM jobid not found in output!")
            job_id = None
    except subprocess.CalledProcessError as e:
        print("[ERROR] Failed to submit SLURM job")
        print(e)
        job_id = None

    return job_id


####################################################################################

# Sample_IDs = ["250206_Occ", "250206_Parie", "250728_Ant", "250728_Post", "250206"]
# TISSUES = ["Tumor", "Tumor", "Tumor", "Tumor", "Blood"]

Sample_IDs = ["SMS_1st", "SMS_2nd", "SMS_3rd", "SMS_4th", "SMS_mets", "SMS"]
TISSUES = ["Tumor", "Tumor", "Tumor", "Tumor", "Tumor", "Blood"]

CONTROL_BAM_PATH = os.path.join(BAM_DIR, "hg38", "Blood", "05.Final_bam", "{}_Blood.bam".format(Sample_IDs[-1]))
CONTROL_SAMPLE = "{}_Blood".format(Sample_IDs[-1])

if os.path.exists(CONTROL_BAM_PATH):     # File이 있어야만 진행
    Sample_ID = Sample_IDs[-1]
    SAMPLE_THRESHOLD = "all"
    DP_THRESHOLD = 10
    ALT_THRESHOLD = 1
    REMOVE_MULTIALLELIC = "True"
    PASS = "True"
    REMOVE_MITOCHONDRIAL_DNA = "True"
    BLACKLIST = "/home/goldpm1/resources/RM+SegDup.bed"
    MIN_BQ = 20
    SEQUENCING = "WGS"

    OUTPUT_VCF_GZ = os.path.join( MUTECT_DIR,  "01.raw", f"{Sample_ID}_multiple.vcf.gz")
    OUTPUT_FMC_PATH = os.path.join(MUTECT_DIR,  "02.PASS", f"{Sample_ID}_multiple.MT2.FMC.vcf")
    OUTPUT_FMC_HF_PATH = os.path.join(MUTECT_DIR,  "02.PASS", f"{Sample_ID}_multiple.MT2.FMC.HF.vcf")
    OUTPUT_FMC_HF_RMBLACK_PATH = os.path.join(MUTECT_DIR,  "02.PASS", f"{Sample_ID}_multiple.MT2.FMC.HF.RMBLACK.vcf")

    folders = [  os.path.dirname(OUTPUT_VCF_GZ), os.path.dirname(OUTPUT_FMC_PATH), os.path.dirname(OUTPUT_FMC_HF_PATH),  os.path.dirname(OUTPUT_FMC_HF_RMBLACK_PATH)   ]
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)


    # 01. Multiple mutect2 call
    gatk_command = "gatk --java-options \"-Xmx48g\" Mutect2 -R {} -normal {} --panel-of-normals {} --germline-resource {} -O {} --tmp-dir {}".\
            format(REF, CONTROL_SAMPLE, PON, gnomad, OUTPUT_VCF_GZ, TMP_PATH)
    for i in range(len(Sample_IDs) - 1):  # 맨 끝은 빼고
        Sample_ID = Sample_IDs[i]
        TISSUE = TISSUES[i]
        gatk_command += " -I {}".format(os.path.join(BAM_DIR, "hg38", TISSUE, "05.Final_bam", f"{Sample_ID}_{TISSUE}.bam"))
    gatk_command += " -I {}".format ( CONTROL_BAM_PATH ) # 마지막에 Blood 추가

    slurm_command = f"sbatch --job-name=MT_11.{Sample_IDs[-1]}_multiple --output={logPath}/11.multipleMT/MT_11.{Sample_IDs[-1]}_multiple.out --error={logPath}/11.multipleMT/MT_11.{Sample_IDs[-1]}_multiple.err \
                                        --cpus-per-task=8 --mem-per-cpu=14GB --time=30:00:00 --partition=cpu --qos nstumor \
                                        --wrap='bash {CURRENT_PATH}/mutect_multiple_pipe_01.call.sh --COMMAND \"{gatk_command}\" --OUTPUT_VCF_GZ {OUTPUT_VCF_GZ}'"
    MT01_JOBID = get_slurm_jobid ( slurm_command )
    print(f"SLURM jobid: {MT01_JOBID}")
    print (slurm_command)



    # # 02. FMC & HF & RMBLACK & VEP (Python version)    
    #--dependency=afterok:{MT01_JOBID}
    slurm_command = f"sbatch  --job-name=MT_02.{Sample_IDs[-1]}_multiple --output={logPath}/02.FMC_HF_RMBLACK/MT_02.{Sample_IDs[-1]}_multiple.out --error={logPath}/02.FMC_HF_RMBLACK/MT_02.{Sample_IDs[-1]}_multiple.err \
                                        --cpus-per-task=5 --mem-per-cpu=14GB --time=6:00:00 --partition=cpu --qos nstumor \
                                        --wrap='bash {CURRENT_PATH}/mutect_pair_pipe_02.FMC_HF_RMBLACK.sh \
                                                    --Sample_ID {Sample_IDs[-1]} --OUTPUT_VCF_GZ {OUTPUT_VCF_GZ} --OUTPUT_FMC_PATH {OUTPUT_FMC_PATH} --OUTPUT_FMC_HF_PATH {OUTPUT_FMC_HF_PATH} \
                                                    --OUTPUT_FMC_HF_RMBLACK_PATH {OUTPUT_FMC_HF_RMBLACK_PATH} \
                                                    --PON {PON} --REF {REF} --gnomad {gnomad} --INTERVAL {INTERVAL} --TMP_PATH {TMP_PATH} \
                                                    --SAMPLE_THRESHOLD {SAMPLE_THRESHOLD} --DP_THRESHOLD {DP_THRESHOLD} --ALT_THRESHOLD {ALT_THRESHOLD} \
                                                    --MIN_BQ {MIN_BQ} --REMOVE_MULTIALLELIC {REMOVE_MULTIALLELIC} --PASS {PASS} --REMOVE_MITOCHONDRIAL_DNA {REMOVE_MITOCHONDRIAL_DNA} --BLACKLIST {BLACKLIST}'"
    MT02_JOBID  = get_slurm_jobid(slurm_command)
    print(f"SLURM jobid: {MT02_JOBID}")
    
    #print (slurm_command)





#qsub_command = f"qsub -pe smp 6 -e {logPath}/11.multipleMT -o {logPath}/11.multipleMT -N MT_11.{Sample_IDs[-1]}_multiple {CURRENT_PATH}/mutect_multiple_pipe_01.call.sh --COMMAND \"{command}\" --OUTPUT_VCF_GZ {OUTPUT_VCF_GZ}"
#print (qsub_command)
#os.system(qsub_command)

# qsub_command = (
#     f'qsub -pe smp 5 '
#     f'-e {logPath}"/02.FMC_HF_RMBLACK" -o {logPath}"/02.FMC_HF_RMBLACK" '
#     f'-hold_jid MT_11.{Sample_IDs[-1]}_multiple '
#     f'-N MT_02.{Sample_ID[-1]}_multiple '
#     f'mutect_pair_pipe_02.FMC_HF_RMBLACK.sh '
#     f'--Sample_ID {Sample_ID[-1]} '
#     f'--OUTPUT_VCF_GZ {OUTPUT_VCF_GZ} '
#     f'--OUTPUT_FMC_PATH {OUTPUT_FMC_PATH} '
#     f'--OUTPUT_FMC_HF_PATH {OUTPUT_FMC_HF_PATH} '
#     f'--OUTPUT_FMC_HF_RMBLACK_PATH {OUTPUT_FMC_HF_RMBLACK_PATH} '
#     f'--PON {PON} --REF {REF} --gnomad {gnomad} --INTERVAL {INTERVAL} --TMP_PATH {TMP_PATH} '
#     f'--SAMPLE_THRESHOLD {SAMPLE_THRESHOLD} --DP_THRESHOLD {DP_THRESHOLD} --ALT_THRESHOLD {ALT_THRESHOLD} '
#     f'--MIN_BQ {MIN_BQ} --REMOVE_MULTIALLELIC {REMOVE_MULTIALLELIC} --PASS {PASS} '
#     f'--REMOVE_MITOCHONDRIAL_DNA {REMOVE_MITOCHONDRIAL_DNA} --BLACKLIST {BLACKLIST}'
# )
# print (qsub_command)
# os.system(qsub_command)