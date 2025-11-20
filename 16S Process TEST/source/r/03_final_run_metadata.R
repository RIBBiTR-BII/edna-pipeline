# Title: Finalize metadata for amplicon sequence processing
# Created by: Cob Staines (cobstainesconsulting@gmail.com)

library(yaml)

# inherit env_config_path
args = commandArgs(trailingOnly = TRUE)
env_config_path = args[1]
#  env_config_path = "runs/2025-11-03_test/output/metadata/config.yml"

# read in config for runDir
config = read_yaml(env_config_path)
# read in run metadata
metadata = read_yaml(paste0(config$run$runDir, "/output/metadata/run_metadata.yml"))
# update end time
end_time = Sys.time()
metadata$post$runtime$endTimestamp = format(end_time, '%Y-%m-%d %H:%M:%S %z')
# calculate duration
start_time = as.POSIXct(metadata$post$runtime$startTimestamp, format = "%Y-%m-%d %H:%M:%S %z")
metadata$post$runtime$duration = paste0(as.character(round(difftime(end_time, start_time, units = "secs"), digits = 0)), " seconds")
# update exit status
metadata$post$exitStatus = dplyr::case_when(
  file.exists(paste0(config$run$runDir, "/output/", config$run$name, "_eDNA_result_tables.xlsx")) ~ "complete",
  dir.exists(paste0(config$run$runDir, "/analysis/s07_classified_taxonomy_vsearch/search_results/")) ~ "incomplete-s07-b",
  dir.exists(paste0(config$run$runDir, "/analysis/s07_classified_taxonomy_cblast/classification/")) ~ "incomplete-s07-a",
  dir.exists(paste0(config$run$runDir, "/analysis/s06_denoised_16S_eDNA/representative_sequences")) ~ "incomplete-s06",
  file.exists(paste0(config$run$runDir, "/analysis/s05_primertrimmed_16S_eDNA_Demux.qzv ")) ~ "incomplete-s05",
  file.exists(paste0(config$run$runDir, "/analysis/s04_leadlag_primertrimmed_16S_eDNA_Demux.qza")) ~ "incomplete-s04",
  file.exists(paste0(config$run$runDir, "/analysis/s03_lead_primertrimmed_16S_eDNA_Demux.qza")) ~ "incomplete-s03",
  file.exists(paste0(config$run$runDir, "/analysis/s02_16S_eDNA_Demux.qzv")) ~ "incomplete-s02",
  file.exists(paste0(config$run$runDir, "/analysis/s01_16S_eDNA_Demux.qza")) ~ "incomplete-s01",
  file.exists(paste0(config$run$runDir, "/output/metadata/run_metadata.yml")) ~ "incomplete-r02",
  file.exists(paste0(config$run$runDir, "/output/metadata/manifest.csv")) ~ "incomplete-r01",
  .default = "incomplete-r00"
)
# write to file
write_yaml(metadata, paste0(config$run$runDir, "/output/metadata/run_metadata.yml"))
