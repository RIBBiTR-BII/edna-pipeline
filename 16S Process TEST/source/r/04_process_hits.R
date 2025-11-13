# Title: Organizing data from NGS for further analysis. 
# Created by: Brandon Hoenig (brandonhoenig@gmail.com)
# Created on: 12 August 2024
# Edited on: 4 November 2024

## Load in Libraries (install them in you need them.)
library(tidyverse)
library(iucnredlist)
library(yaml)

# 
# # # manual runs
setwd("16S Process TEST")
env_config_path = "runs/2025-11-07_panama/output/metadata/config.yml"

# read in config file
config = read_yaml(env_config_path)

hit_table_raw = read_csv(paste0(config$run$runDir, "/output/", config$run$name, "_eDNA_hits_joined.csv"))

hit_table_clean = hit_table_raw %>%
  filter(!is.na(class)) %>%
  mutate(
    species_simple = gsub(" .*", "", species),
    species_simple = if_else(species_simple %in% c("sp.",
                                                   "aff.",
                                                   "cf."),
                             NA,
                             species_simple),
    scientific_name = case_when(
      !is.na(genus) & !is.na(species_simple) ~ paste0(genus, " ", species_simple),
      !is.na(genus) & is.na(species_simple) ~ genus,
      .default = NA)
  )

species_unique = hit_table_clean %>%
  select(genus, species_simple, scientific_name) %>%
  filter(!is.na(species_simple)) %>%
  distinct() %>%
  arrange(scientific_name) %>%
  mutate()

# IUCN search for thos in Panama

api = init_api(Sys.getenv("iucn_token"))

# init empty df
abn = assessments_by_name(api, genus = species_unique$genus[1], species = species_unique$species_simple[1]) %>%
  mutate(genus = species_unique$genus[1],
         species = species_unique$species_simple[1]) %>%
  filter(FALSE) %>%
  select(genus,
         species,
         everything())

for (ii in 112:nrow(species_unique)) {
  print(ii)
  ii_genus = species_unique$genus[ii]
  ii_species = species_unique$species_simple[ii]
  
  # pull data
  new_row <- tryCatch({
    # Your main code block
    assessments_by_name(api, genus = ii_genus, species = ii_species) %>%
      mutate(genus = ii_genus,
             species = ii_species) %>%
      arrange(desc(year_published)) %>%
      slice_head(n = 1)
  }, error = function(e) {
    message("Error occurred, skipping this row.")
    NULL
  })
  
  
  # bind
  if (!is.null(new_row) && nrow(new_row) > 0) {
    abn <- bind_rows(abn, new_row)
  }
  # sleep
  Sys.sleep(0.5)
}

write_csv(abn, paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments.csv"))
# abn = read_csv(paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments.csv"))


a_data = assessment_data_many(api,
                              abn$assessment_id,
                              wait_time = 0.5)
saveRDS(a_data, file = paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments_pulled.RData"))
# load(paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments_pulled.RData"))

assessment = a_data[[285]]

found_in_panama = function(assessment) {
  locations = assessment$locations %>%
    filter(code == "PA")
    
}

assessment$assessment_id$value

loc_panama = map_dfr(a_data, function(x) {
  x$location %>%
    filter(code == "PA") %>%
    mutate(assessment_id = x$assessment_id$value,
           seasonality = as.character(seasonality),
           formerlyBred = as.character(formerlyBred))
}) %>%
  mutate(found_in_panama = TRUE)

hit_panama = hit_table_clean %>%
  left_join(abn, by = c("genus", "species_simple" = "species")) %>%
  left_join(loc_panama, by = "assessment_id") %>%
  mutate(found_in_panama = if_else(is.na(found_in_panama) & !is.na(assessment_id),
                                   FALSE,
                                   found_in_panama)) %>%
  filter(found_in_panama | is.na(found_in_panama)) %>%
  select(any_of(colnames(hit_table_clean)),
         found_in_panama)

hit_panama_consensus_value = hit_panama %>%
  filter(percent_identical >= 0.75) %>%
  group_by(asv_id, method) %>%
  mutate(hits = n()) %>%
  ungroup() %>%
  group_by(asv_id, method, class, order, family, genus, species) %>%
  mutate(taxon_hits = n()) %>%
  ungroup() %>%
  mutate(consensus_value = taxon_hits / hits) %>%
  arrange(desc(asv_total_count),
          asv_id,
          method,
          desc(consensus_value),
          percent_identical,
          seq_overlap) %>%
  group_by(asv_id, method, class, order, family, genus, species) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(desc(asv_total_count),
          asv_id,
          method)

hit_panama_vsearch_consensus = hit_panama_consensus %>%
  filter(method == "vsearch",
         consensus_value >= 0.5)

# consensus blast results
consensus_table_blast =
  read.table(paste0(config$run$runDir, "/analysis/s07_classified_taxonomy_blast/classification/", 
                    list.files(paste0(config$run$runDir, '/analysis/s07_classified_taxonomy_blast/classification')), 
                    "/data/taxonomy.tsv"),
             header = TRUE,
             sep = '\t') %>%
  separate(Taxon, into = c("k", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate(method = "blast",
         across(everything(), ~ str_remove(., "^[a-z]__"))) %>%
  rename(asv_id = Feature.ID,
         kingdom = k,
         phylum = p,
         class = c,
         order = o,
         family = f,
         genus = g,
         species = s,
         consensus = Consensus) %>%
  select(-kingdom,
         -phylum)

# decision tree
accept_vsearch = hit_table_raw %>%
  select(asv_id,
         asv_total_count,
         sample_count,
         samples,
         sequence) %>%
  distinct() %>%
  left_join(hit_panama_vsearch_consensus %>%
              select(-asv_total_count,
                     -sample_count,
                     -samples,
                     -sequence), by = "asv_id") %>%
  select(-consensus) %>%
  rename(consensus = consensus_value)

vsearch_assigned = sum(!is.na(accept_vsearch$class)) / nrow(accept_vsearch)

accept_cblast = hit_table_raw %>%
  select(asv_id,
         asv_total_count,
         sample_count,
         samples,
         sequence) %>%
  distinct() %>%
  left_join(consensus_table_blast, by = "asv_id")

cblast_assigned = sum(!is.na(accept_cblast$class)) / nrow(accept_cblast)

discrepancies = accept_vsearch %>%
  select(any_of(colnames(accept_cblast))) %>%
  left_join(accept_cblast %>%
              select(-asv_total_count,
                     -sample_count,
                     -samples,
                     -sequence), by = "asv_id") %>%
  filter((species.x != species.y) | (genus.x != genus.y)) 

unassigned = hit_panama_consensus_value %>%
  anti_join(accept_vsearch %>%
              filter(!is.na(class)), by = "asv_id")



# limit to amphibean hits
