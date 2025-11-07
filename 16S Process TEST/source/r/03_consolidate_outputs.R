# Title: Organizing data from NGS for further analysis. 
# Created by: Brandon Hoenig (brandonhoenig@gmail.com)
# Created on: 12 August 2024
# Edited on: 4 November 2024

## Load in Libraries (install them in you need them.)
library(tidyverse)
library(phylotools)
library(readxl)
library(openxlsx)
library(yaml)

# # manual runs
# env_config_path = "runs/2025-11-06_penn/metadata/config.yml"

# inherit env_config_path
args = commandArgs(trailingOnly = TRUE)
env_config_path = args[1]

# read in config file
config = read_yaml(env_config_path)

## Read in the .biom feature table. 
feature_table =
  read_tsv(paste0(config$run$runDir, "/analysis/s05_denoised_16S_eDNA/sample_table/feature-table.tsv"), skip = 1) %>%
  as.data.frame() %>%
  rename("asv_id" = "#OTU ID") %>%
  arrange(asv_id)

feature_long = feature_table %>%
  pivot_longer(cols = -asv_id,
               names_to = "sample_id",
               values_to = "asv_count") %>%
  filter(asv_count != 0) %>%
  arrange(sample_id) %>%
  select(sample_id,
         asv_id,
         asv_count)

asv_total_count = feature_long %>%
  group_by(asv_id) %>%
  summarise(asv_total_count = sum(asv_count)) %>%
  arrange(desc(asv_total_count))

## Read in taxonomic classifications (vsearch global results)
classification_table =
  read.table(paste0(config$run$runDir, "/analysis/s06_classified_taxonomy_vsearch/search_results/", 
                    list.files(paste0(config$run$runDir, '/analysis/s06_classified_taxonomy_vsearch/search_results')), 
                    "/data/blast6.tsv"),
             header = FALSE, 
             sep = '\t') %>%
  rename(asv_id = V1,
         taxon_id = V2,
         percent_identical = V3,
         seq_overlap = V4,
         seq_mismatch = V5,
         gapopen_count = V6,
         q_start = V7,
         q_end = V8,
         s_start = V9,
         s_end = V10,
         e_value = V11,
         bitscore = V12) %>%
  arrange(asv_id)

# Read in taxonomies and filter to flagged taxa
taxonomy_table =
  read.table(paste0(config$taxonomy$classifierDir, "/Vertebrata16S_derep1_taxa_extracted/", 
                    list.files(paste0(config$taxonomy$classifierDir, '/Vertebrata16S_derep1_taxa_extracted')), 
                    "/data/taxonomy.tsv"),
             header = TRUE, 
             sep = '\t',
             quote = '') %>%
  rename(taxon_id = Feature.ID,
         taxon = Taxon) %>%
  filter(taxon_id %in% unique(classification_table$taxon_id)) %>%
  separate(taxon, into = c("k", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate(across(everything(), ~ str_remove(., "^[a-z]__"))) %>%
  arrange(taxon_id)

# Read in DNA sequences 
sequence_table =
  read.fasta(file = paste0(config$run$runDir, "/analysis/s05_denoised_16S_eDNA/representative_sequences/", 
                           list.files(paste0(config$run$runDir, '/analysis/s05_denoised_16S_eDNA/representative_sequences')), 
                           "/data/dna-sequences.fasta"))%>%
  rename(asv_id = seq.name, 
         sequence = seq.text) %>%
  arrange(asv_id)

wide_table = asv_total_count %>%
  left_join(sequence_table, by = "asv_id") %>%
  left_join(classification_table, by = "asv_id") %>%
  left_join(taxonomy_table, by = "taxon_id")

write.csv(wide_table, paste0(config$run$runDir, "/output/", config$run$name, "_eDNA_hits_joined.csv"))

## Write each table for documentation.  Proceed with analysis in R or in Excel.
list_of_datasets <- list("Feature Table" = feature_table, 
                         "Taxonomy Table" = taxonomy_table,
                         "Classification Table" = classification_table,
                         "Sequence Table" = sequence_table)

write.xlsx(list_of_datasets, file = paste0(config$run$runDir, "/output/", config$run$name, "_eDNA_result_tables.xlsx"))
