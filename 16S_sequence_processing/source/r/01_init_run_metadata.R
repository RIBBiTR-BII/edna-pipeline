# Title: Initiate metadata for amplicon sequence processing
# Created by: Cob Staines (cobstainesconsulting@gmail.com)

library(yaml)

# # manual runs
# env_config_path = "runs/test_run_01/output/metadata/config.yml"
# setwd("16S_sequence_processing")

# inherit env_config_path
args = commandArgs(trailingOnly = TRUE)
env_config_path = args[1]

metadata = list()
# read in config file, use as basis for metadata
metadata$config = read_yaml(env_config_path)

# check if classifier specified
if (!is.null(metadata$config$taxonomy$classifierDir)) {
  # scrape classifier metadata
  classifier_md = read_yaml(paste0(metadata$config$taxonomy$classifierDir, "/classifier_metadata.yml"))
  metadata$post$classifier = classifier_md$classifier
}

# scrape manifest file
manifest = read.csv(paste0(metadata$config$run$runDir, "/output/metadata/manifest.csv"))
metadata$post$sequences$sampleCount = length(unique(manifest$sample.id))
metadata$post$sequences$sampleIDFirst = head(manifest$sample.id, 1)
metadata$post$sequences$sampleIDLast = tail(manifest$sample.id, 1)
# add start and end time
metadata$post$runtime$startTimestamp = format(Sys.time(), '%Y-%m-%d %H:%M:%S %z')
metadata$post$runtime$endTimestamp = NA
# add exit status
metadata$post$exitStatus = "incomplete"
#write to file
write_yaml(metadata, paste0(metadata$config$run$runDir, "/output/metadata/run_metadata.yml"))
