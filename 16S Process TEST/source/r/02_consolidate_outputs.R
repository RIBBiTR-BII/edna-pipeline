# Title: Organizing data from NGS for further analysis. 
# Created by: Brandon Hoenig (brandonhoenig@gmail.com) & Cob Staines (cobstainesconsulting@gmail.com)

## Load in Libraries (install as needed)
library(tidyverse)
library(phylotools)
library(readxl)
library(openxlsx)
library(yaml)

# 
# # # manual runs
# env_config_path = "runs/methods_2025-12-29/output/metadata/config.yml"
# setwd("16S Process TEST")

# inherit env_config_path
args = commandArgs(trailingOnly = TRUE)
env_config_path = args[1]

# read in config file
config = read_yaml(env_config_path)

## Read in the .biom feature table. 
feature_table =
  read_tsv(paste0(config$run$runDir, "/analysis/s06_denoised_16S_eDNA/sample_table/feature-table.tsv"), skip = 1) %>%
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

asv_totals = feature_long %>%
  group_by(asv_id) %>%
  summarise(asv_total_count = sum(asv_count),
            sample_count = n_distinct(sample_id),
            samples = paste(sort(unique(sample_id)), collapse = ", ")) %>%
  arrange(desc(asv_total_count))

## Read in taxonomic classifications (vsearch global results)
hit_table_vsearch =
  read.table(paste0(config$run$runDir, "/analysis/s07_classified_taxonomy_vsearch/search_results/", 
                    list.files(paste0(config$run$runDir, '/analysis/s07_classified_taxonomy_vsearch/search_results')), 
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
  arrange(asv_id) %>%
  mutate(method = "vsearch")

hit_table_blast =
  read.table(paste0(config$run$runDir, "/analysis/s07_classified_taxonomy_blast/search_results/", 
                    list.files(paste0(config$run$runDir, '/analysis/s07_classified_taxonomy_blast/search_results')), 
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
  arrange(asv_id) %>%
  mutate(method = "blast")

consensus_table_blast =
  read.table(paste0(config$run$runDir, "/analysis/s07_classified_taxonomy_blast/classification/", 
                                          list.files(paste0(config$run$runDir, '/analysis/s07_classified_taxonomy_blast/classification')), 
                                          "/data/taxonomy.tsv"),
             header = TRUE,
             sep = '\t') %>%
  separate(Taxon, into = c("k", "p", "c", "o", "f", "g", "s"), sep = ";", fill = "right") %>%
  mutate(method = "blast",
         across(everything(), ~ str_remove(., "^[a-z]__")),
         k = if_else(k == "Unassigned", NA, k)) %>%
  rename(asv_id = Feature.ID,
         kingdom = k,
         phylum = p,
         class = c,
         order = o,
         family = f,
         genus = g,
         species = s,
         consensus = Consensus)

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
  filter(taxon_id %in% unique(hit_table_vsearch$taxon_id) |
           taxon_id %in% unique(hit_table_blast$taxon_id)) %>%
  separate(taxon, into = c("k", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate(across(everything(), ~ str_remove(., "^[a-z]__"))) %>%
  rename(kingdom = k,
         phylum = p,
         class = c,
         order = o,
         family = f,
         genus = g,
         species = s) %>%
  arrange(taxon_id)

# Read in DNA sequences 
sequence_table =
  read.fasta(file = paste0(config$run$runDir, "/analysis/s06_denoised_16S_eDNA/representative_sequences/", 
                           list.files(paste0(config$run$runDir, '/analysis/s06_denoised_16S_eDNA/representative_sequences')), 
                           "/data/dna-sequences.fasta"))%>%
  rename(asv_id = seq.name, 
         sequence = seq.text) %>%
  arrange(asv_id)

wide_table = asv_totals %>%
  left_join(sequence_table, by = "asv_id") %>%
  left_join(bind_rows(hit_table_blast,
                      hit_table_vsearch), by = "asv_id") %>%
  left_join(taxonomy_table, by = "taxon_id") %>%
  left_join(consensus_table_blast, by = c("asv_id", "method", "kingdom", "phylum", "class", "order", "family", "genus", "species")) %>%
  select(any_of(colnames(asv_totals)),
         any_of(colnames(sequence_table)),
         method,
         consensus,
         any_of(colnames(taxonomy_table)),
         any_of(colnames(hit_table_vsearch)),
         -kingdom,
         -phylum,
         -e_value,
         -bitscore) %>%
  arrange(desc(asv_total_count), asv_id, method)

write.csv(wide_table, paste0(config$run$runDir, "/output/", config$run$name, "_eDNA_hits_joined.csv"), row.names = FALSE)
write.csv(taxonomy_table, paste0(config$run$runDir, "/output/", config$run$name, "_taxonomy_table_filtered.csv"), row.names = FALSE)

## Write each table for documentation.  Proceed with analysis in R or in Excel.
list_of_datasets <- list("Feature Table" = feature_table, 
                         "Taxonomy Table" = taxonomy_table,
                         "Sequence Table" = sequence_table,
                         "Hits vsearch" = hit_table_vsearch,
                         "Hits BLAST" = hit_table_blast,
                         "Consensus BLAST" = consensus_table_blast)

write.xlsx(list_of_datasets, file = paste0(config$run$runDir, "/output/", config$run$name, "_eDNA_result_tables.xlsx"))
