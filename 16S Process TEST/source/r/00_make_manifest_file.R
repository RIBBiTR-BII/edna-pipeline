library(tidyverse)
library(yaml)

# inherit env_config_path
args = commandArgs(trailingOnly = TRUE)
env_config_path = args[1]

# read in config file
config = read_yaml(env_config_path)

# create manifest file in tabular form
tibble(`absolute-filepath` = c(list.files(paste0(config$run$runDir, "/sequences"), recursive = TRUE))) %>%
  filter(str_detect(`absolute-filepath`, "\\.fastq(\\.gz)?$")) %>%
  mutate(`seq_id` = str_match(basename(`absolute-filepath`), "([A-Z]{2}[0-9]{3})-[0-9A-Za-z-]+_S[0-9]+_L[0-9]+_R[12]")[,2],
         `sample-id` = str_match(basename(`absolute-filepath`), "-([A-Za-z]+-[0-9A-Z-]*)_S[0-9]+_L[0-9]+_R[12]")[,2],
         `absolute-filepath` = paste0("$PWD/sequences/", `absolute-filepath`),
         `direction` = if_else(str_detect(`absolute-filepath`, "_R1_"), "forward", 'reverse'),
         `sample-id` = if_else(is.na(`sample-id`), `seq_id`, `sample-id`),
         `sample-id` = str_replace(`sample-id`, "-", "_")) %>%
  select(`sample-id`, `absolute-filepath`, `direction`) %>%
  # write as a csv and place in metadata folder
  write_csv(., paste0(config$run$runDir, "/metadata/manifest.csv"))
