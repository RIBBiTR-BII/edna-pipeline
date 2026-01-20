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
env_config_path = "runs/methods_2025-12-29/output/metadata/config.yml"

# read in config file
config = read_yaml(env_config_path)

data_filtered = read.csv(paste0(config$run$runDir, "/output/", config$run$name, "_all_samples_nc_filter.csv"))

data_amphibia_species = data_filtered %>%
  filter(class == "Amphibia",
         taxonomic_level == "species")

sample_species_asv = data_amphibia_species %>%
  group_by(site, genus, species) %>%
  mutate(pop_haplotype_n = n(),
         pop_asv_total = sum(asv_count),
         pop_haplotype_frequency = asv_count / pop_asv_total) %>%
  ungroup() %>%
  group_by(site, date, genus, species) %>%
  mutate(sample_haplotype_n = n(),
         sample_asv_total = sum(asv_count),
         sample_haplotype_frequency = asv_count / sample_asv_total) %>%
  ungroup() %>%
  arrange(desc(sample_haplotype_n),
          site,
          date)

estimate_species = sample_species_asv %>%
  group_by(site, date, genus, species) %>%
  summarize(pop_haplotype_n = first(pop_haplotype_n),
            sample_n = n_distinct(illumina_sample_name),
            sample_asv_total = first(sample_asv_total),
            sample_haplotype_n = n(),
            sample_asv_total = sum(asv_count),
            big_n = n()/sum((sample_haplotype_frequency - pop_haplotype_frequency)^2/(pop_haplotype_frequency*(1 - pop_haplotype_frequency))),
            .groups = "drop") %>%
  arrange(site, species, date)
