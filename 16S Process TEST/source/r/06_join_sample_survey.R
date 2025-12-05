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
db_env = tbl(dbcon, Id("survey_data", "environmental"))
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


all_ids = feature_long = feature_table %>%
  pivot_longer(cols = -asv_id,
               names_to = "illumina_sample_raw",
               values_to = "asv_count") %>%
  mutate(illumina_sample_name = gsub("-B$", "", illumina_sample_raw),
         illumina_sample_name = gsub("(_\\d{2}-)(\\d{2}$)", "\\10\\2", illumina_sample_name)) %>%
  select(illumina_sample_raw,
         illumina_sample_name) %>%
  distinct()

feature_long = feature_table %>%
  pivot_longer(cols = -asv_id,
               names_to = "illumina_sample_raw",
               values_to = "asv_count") %>%
  filter(asv_count != 0) %>%
  mutate(illumina_sample_name = gsub("-B$", "", illumina_sample_raw),
         illumina_sample_name = gsub("(_\\d{2}-)(\\d{2}$)", "\\10\\2", illumina_sample_name)) %>%
  arrange(illumina_sample_name) %>%
  select(illumina_sample_name,
         illumina_sample_raw,
         asv_id,
         asv_count) %>%
  full_join(all_ids, by = c("illumina_sample_raw", "illumina_sample_name"))

unique_samples = feature_long %>%
  select(illumina_sample_name) %>%
  distinct() %>%
  arrange(illumina_sample_name)

# panama survey data
edna_panama = db_edna %>%
  left_join(db_survey, by = "survey_id") %>%
  left_join(db_visit, by = "visit_id") %>%
  left_join(db_site, by = "site_id") %>%
  left_join(db_region, by = "region_id") %>%
  left_join(db_country, by = "country_id") %>%
  filter(country == "panama") %>%
  collect()

samps_joined = unique_samples %>%
  left_join(edna_panama, by = "illumina_sample_name")

env_panama = db_env %>%
  left_join(db_survey, by = "survey_id") %>%
  left_join(db_visit, by = "visit_id") %>%
  filter(visit_id %in% samps_joined$visit_id) %>%
  collect() %>%
  arrange(water_time) %>%
  group_by(visit_id, transect) %>%
  slice_head(n = 1) %>%
  ungroup()

samps_env_joined = samps_joined %>%
  left_join(env_panama %>%
              select(visit_id,
                     transect,
                     any_of(colnames(db_env))), by = c("visit_id", "transect")) %>%
  select(illumina_sample_name,
         site,
         transect,
         date,
         collection_time,
         collection_location,
         filter_date,
         filter_method,
         filter_size_um,
         filter_volume_ml,
         filter_start_time,
         filter_end_time,
         negative_control,
         negative_control_group_id,
         comments_edna,
         comments_lab,
         edna_transect_m,
         edna_latitude,
         edna_longitude,
         air_time,
         air_temp_c,
         wind_speed_m_s,
         cloud_cover_percent,
         water_time,
         water_temp_c,
         p_h,
         tds_ppm,
         salinity_ppt,
         conductivity_us_cm,
         visit_id) %>%
  arrange(site, transect, date, collection_location, filter_method, filter_size_um)

write.csv(feature_long, paste0(config$run$runDir, "/output/", config$run$name, "_feature_table_long.csv"), row.names = FALSE)

write.csv(samps_env_joined, paste0(config$run$runDir, "/output/", config$run$name, "_samples_environment.csv"), row.names = FALSE)

# complete output
sample_output = feature_long %>%
  left_join(asv_classified, by = "asv_id") %>%
  left_join(samps_env_joined, by = "illumina_sample_name") %>%
  arrange(site, transect, date, collection_location, filter_method, filter_size_um, desc(asv_count))

write.csv(sample_output, paste0(config$run$runDir, "/output/", config$run$name, "_samples_environment_asv_taxonomy.csv"), row.names = FALSE)

# drop taxa in negative controls

## within negative control group, list all identified taxa in negative control
nc_filter_group = sample_output %>%
  filter(negative_control) %>%
  filter(!is.na(method))

nc_filter_global = sample_output %>%
  filter(grepl("^[A-Z]{2}\\d{3}$", illumina_sample_name)) %>%
  filter(!is.na(method))

# antijoin with each group at highest taxon level
sample_nc_filter = sample_output %>%
  anti_join(nc_filter_group, by = c("asv_id", "negative_control_group_id")) %>%
  anti_join(nc_filter_global, by = c("asv_id")) %>%
  filter(!is.na(method))

sample_filter_all = sample_nc_filter %>%
  bind_rows(samps_env_joined %>%
              filter(!is.na(site),
                     !negative_control)) %>%
  mutate(asv_count = if_else(is.na(asv_count), 0, asv_count))
  
write.csv(sample_filter_all, paste0(config$run$runDir, "/output/", config$run$name, "_all_samples_nc_filter.csv"), row.names = FALSE)
