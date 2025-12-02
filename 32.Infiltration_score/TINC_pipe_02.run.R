library ( argparse )
# "TINC"  추가
tinc_lib_path <- "/home/goldpm1/miniforge3/envs/TINC/lib/R/library"
.libPaths(c(tinc_lib_path, .libPaths()))
library ( ggplot2 )
library ( crayon )
library ( TINC )
library ( mobster )

# /home/goldpm1/tools/TINC/TINC/R/analyses.R에서   pi_cutoff = 0.02,  N_cutoff = 10 -> 2으로 변경한 후 저장하고
#  R CMD INSTALL /home/goldpm1/tools/TINC/TINC  --library="/home/goldpm1/miniforge3/envs/TINC/lib/R/library" 로 재설치하자



conda_env <- Sys.getenv("CONDA_DEFAULT_ENV")
if (conda_env != "TINC") {
  warning(paste0("This script should be run in the 'TINC' conda environment (currently: '", conda_env, "'). Activate with: conda activate TINC"))
} else {
  paste0("TINC conda environment is activated")
}


# Input by argparser
parser <- ArgumentParser()
parser$add_argument("--Sample_ID", default = "230419")
#parser$add_argument("--Sample_ID", default = "230419")
parser$add_argument("--MIN_VAF", default = 0.03 )
parser$add_argument("--N_cutoff", default = 4 )
parser$add_argument("--pi_cutoff", default = 0.03)
parser$add_argument("--SNV_INPUT", default = "/data/project/Meningioma/32.Infiltration_score/01.TINC/01.Input/02.SNV/230419.tsv")
parser$add_argument("--CNV_INPUT", default = "/data/project/Meningioma/32.Infiltration_score/01.TINC/01.Input/01.CNV/230419.tsv")
parser$add_argument("--OUTPUT_PATH", default = "/data/project/Meningioma/32.Infiltration_score/01.TINC/02.Output/230419/TINC_output.tsv")
parser$add_argument("--OUTPUT_PLOT_DIR", default = "/data/project/Meningioma/32.Infiltration_score/01.TINC/02.Output/230419")
args <- parser$parse_args()

Sample_ID = args$Sample_ID
SNV_INPUT =  args$SNV_INPUT
CNV_INPUT = args$CNV_INPUT
MIN_VAF = as.numeric(args$MIN_VAF)
N_CUTOFF = args$N_cutoff
PI_CUTOFF = args$pi_cutoff
OUTPUT_PATH = args$OUTPUT_PATH
OUTPUT_PLOT_DIR = args$OUTPUT_PLOT_DIR

parent_dir <- dirname(OUTPUT_PATH)
if (!file.exists(parent_dir)) {
  dir.create(parent_dir, recursive = TRUE)
}

args=commandArgs(trailingOnly=TRUE)



somatic_data = read.table(SNV_INPUT, header = TRUE, sep = "\t")

if ( file.exists(CNV_INPUT) == TRUE) {
    cnv_data = read.table(CNV_INPUT, header = TRUE, sep = "\t")
    TINC_fit = autofit (  input = somatic_data,  cna   = cnv_data,  VAF_range_tumour = c(MIN_VAF, 0.7),  cutoff_miscalled_clonal = 0.7, N_cutoff = N_CUTOFF, pi_cutoff = PI_CUTOFF,  FAST  = TRUE )             # FAST: parallel processing
} else{
    paste0("Cannot find CNV_INPUT file: ", CNV_INPUT)
    #   ,  cutoff_miscalled_clonal = 0.6,
    TINC_fit = autofit (  input = somatic_data,  cna = NULL,   cutoff_miscalled_clonal = 0.7, FAST  = TRUE )             # FAST: parallel processing
}




out <- TINC::to_string(TINC_fit)
if (!is.null(out)) {
  results_df = data.frame( sample = Sample_ID,TumorInNormal = out$TIN,  TumorInTumor = out$TIT,  stringsAsFactors = FALSE )
  write.table(results_df, file = OUTPUT_PATH, sep = "\t", row.names = FALSE, quote = FALSE)
}




print ( TINC_fit )



# Save TINC plot as a PNG file using ggplot2 and ggsave
plot_obj <- plot(TINC_fit)
ggsave(  file.path( OUTPUT_PLOT_DIR, "TINC_plot.jpg" )  , plot = plot_obj, width = 8, height = 6, dpi = 300)

plot_full_page_report <- TINC::plot_full_page_report(TINC_fit)
ggsave(  file.path( OUTPUT_PLOT_DIR, "TINC_full_page_report.jpg" )  , plot = plot_full_page_report, width = 16, height = 12, dpi = 300)
