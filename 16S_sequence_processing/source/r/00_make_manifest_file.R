# Title: Make manifest file
# Created by: Cob Staines (cobstainesconsulting@gmail.com)

library(tidyverse)
library(yaml)

# # manual runs
# env_config_path = "runs/test_run_01/output/metadata/config.yml"
# setwd("16S_sequence_processing")

# inherit env_config_path
args = commandArgs(trailingOnly = TRUE)
env_config_path = args[1]
qc_run = args[2]

qc_run_lgl <- !is.na(qc_run) && isTRUE(as.logical(qc_run))

if (qc_run_lgl) {
  sequence_dir = dirname(env_config_path)
  out_dir = paste0(env_config_path, "/manifest.csv")
  filepath_prefix = "$PWD/"
} else {
  config = read_yaml(env_config_path)
  sequence_dir = paste0(config$run$runDir, "/sequences")
  out_dir = paste0(config$run$runDir, "/output/metadata/manifest.csv")
  filepath_prefix = "$PWD/sequences/"
}

# create manifest file in tabular form
peace = tibble(`absolute-filepath` = c(list.files(sequence_dir, recursive = TRUE))) %>%
  filter(str_detect(`absolute-filepath`, "\\.fastq(\\.gz)?$")) %>%
  mutate(`seq_id` = str_match(basename(`absolute-filepath`), "([A-Z]{2}[0-9]{3})-[0-9A-Za-z-]+_S[0-9]+_L[0-9]+_R[12]")[,2],
         `sample-id` = str_match(basename(`absolute-filepath`), "-([A-Za-z]+-[0-9A-Za-z-]*)_S[0-9]+_L[0-9]+_R[12]")[,2],
         `absolute-filepath` = paste0(filepath_prefix, `absolute-filepath`),
         `direction` = if_else(str_detect(`absolute-filepath`, "_R1_"), "forward", 'reverse'),
         `sample-id` = if_else(is.na(`sample-id`), `seq_id`, `sample-id`),
         `sample-id` = str_replace(`sample-id`, "-", "_")) %>%
  select(`sample-id`, `absolute-filepath`, `direction`) %>%
  # write as a csv and place in metadata folder
  write_csv(., out_dir)

