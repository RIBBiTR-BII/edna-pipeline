# Title: Process taxonomic hits to identify species likely found in Panama
# Created by: Cob Staines (cobstainesconsulting@gmail.com)

## Load in Libraries (install them in you need them.)
library(tidyverse)
library(yaml)
library(dbplyr)
library(RPostgres)
library(DBI)
library(ribbitrrr)


# manual runs
# setwd("16S Process TEST")
env_config_path = "runs/2025-11-07_panama/output/metadata/config.yml"

# read in config file
config = read_yaml(env_config_path)

# establish database connection
dbcon = hopToDB("ribbitr")

# db table pointers
db_edna = tbl(dbcon, Id("survey_data", "edna"))
db_survey = tbl(dbcon, Id("survey_data", "survey"))
db_visit = tbl(dbcon, Id("survey_data", "visit"))
db_site = tbl(dbcon, Id("survey_data", "site"))
db_region = tbl(dbcon, Id("survey_data", "region"))
db_country = tbl(dbcon, Id("survey_data", "country"))

# classified ASVs
asv_classified = read_csv(paste0(config$run$runDir, "/output/", config$run$name, "_asv_classified.csv"))

# samples
feature_table =
  read_tsv(paste0(config$run$runDir, "/analysis/s06_denoised_16S_eDNA/sample_table/feature-table.tsv"), skip = 1) %>%
  as.data.frame() %>%
  rename("asv_id" = "#OTU ID") %>%
  arrange(asv_id)

feature_long = feature_table %>%
  pivot_longer(cols = -asv_id,
               names_to = "illumina_sample_name",
               values_to = "asv_count") %>%
  filter(asv_count != 0) %>%
  arrange(illumina_sample_name) %>%
  select(illumina_sample_name,
         asv_id,
         asv_count)

unique_samples = feature_long %>%
  select(illumina_sample_name) %>%
  distinct()

# panama survey data
edna_panama = db_edna %>%
  left_join(db_survey, by = "survey_id") %>%
  left_join(db_visit, by = "visit_id") %>%
  left_join(db_site, by = "site_id") %>%
  left_join(db_region, by = "region_id") %>%
  left_join(db_country, by = "country_id") %>%
  filter(country == "panama") %>%
  select(any_of(colnames(db_edna)),
         transect,
         site,
         date) %>%
  collect() %>%
  # fix sample name format
  mutate(illumina_sample_name = gsub("edna_", "", illumina_sample_name, ignore.case = TRUE),
         illumina_sample_name = gsub("(PAN_\\d{2})_", "\\1-", illumina_sample_name))

samples_joined = unique_samples %>%
  left_join(edna_panama, by = "illumina_sample_name")
