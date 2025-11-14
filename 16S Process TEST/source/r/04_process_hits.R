# Title: Organizing data from NGS for further analysis. 
# Created by: Brandon Hoenig (brandonhoenig@gmail.com)
# Created on: 12 August 2024
# Edited on: 4 November 2024

## Load in Libraries (install them in you need them.)
library(tidyverse)
library(iucnredlist)
library(rgbif)
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
      .default = NA),
    p_overlap = seq_overlap / s_end,
    p_identical = percent_identical / 100
  )

species_unique = hit_table_clean %>%
  select(genus, species_simple, scientific_name) %>%
  filter(!is.na(species_simple)) %>%
  distinct() %>%
  arrange(scientific_name) %>%
  mutate()

# pull from gbif
gbif_backbone = name_backbone_checklist(species_unique %>%
                                  rename(scientificName = scientific_name,
                                         species = species_simple) %>%
                                  mutate(kingsom = "Animalia"))

rownames(gbif_backbone) = NULL
rownames(species_unique) = NULL

gbif_backbone_key = bind_cols(species_unique %>%
                                rename(species_q = species_simple,
                                       genus_q = genus,
                                       scientific_name_q = scientific_name),
                              gbif_backbone) %>%
  mutate(panama_occ_count = NA)

print("Searching GBIF for species occurences in Panama")
pb = txtProgressBar(min = 1, max = nrow(gbif_backbone_key), style = 3)
for (ii in 1:nrow(gbif_backbone_key)) {
  if (!is.na(gbif_backbone_key$usageKey[ii])) {
    gbif_backbone_key$panama_occ_count[ii] = occ_count(taxonKey = gbif_backbone_key$usageKey[ii],
                                                   country = "PA")
    Sys.sleep(0.1)
    setTxtProgressBar(pb, ii)
  }
}
close(pb)

write_csv(gbif_backbone_key, paste0(config$run$runDir, "/output/", config$run$name, "_GBIF_panama_occurence.csv"))
# gbif_backbone_key = read_csv(paste0(config$run$runDir, "/output/", config$run$name, "_GBIF_panama_occurence.csv"))

gbif_scarce = gbif_backbone_key %>%
  filter(panama_occ_count != 0,
         panama_occ_count < 20,
         !is.na(panama_occ_count))

# cross-reference with IUCN
api = init_api(Sys.getenv("iucn_token"))

# init empty df
abn = assessments_by_name(api, genus = gbif_scarce$genus_q[1], species = gbif_scarce$species_q[1]) %>%
  mutate(genus_q = gbif_scarce$genus_q[1],
         species_q = gbif_scarce$species_q[1]) %>%
  filter(FALSE) %>%
  select(genus_q,
         species_q,
         everything())

print("Cross-referencing with IUCN for scarce species")
pb = txtProgressBar(min = 1, max = nrow(gbif_scarce), style = 3)
for (ii in 1:nrow(gbif_scarce)) {
  ii_genus = gbif_scarce$genus_q[ii]
  ii_species = gbif_scarce$species_q[ii]
  setTxtProgressBar(pb, ii)
  
  # pull data
  new_row <- tryCatch({
    # Your main code block
    assessments_by_name(api, genus = ii_genus, species = ii_species) %>%
      mutate(genus_q = ii_genus,
             species_q = ii_species) %>%
      arrange(desc(year_published)) %>%
      slice_head(n = 1)
  }, error = function(e) {
    # fail silently
    NULL
  })
  
  
  # bind
  if (!is.null(new_row) && nrow(new_row) > 0) {
    abn = bind_rows(abn, new_row)
  }
  # sleep
  Sys.sleep(0.5)
}
close(pb)

write_csv(abn, paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments.csv"))
# abn = read_csv(paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments.csv"))


a_data = assessment_data_many(api,
                              abn$assessment_id,
                              wait_time = 0.5)
saveRDS(a_data, file = paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments_pulled.RData"))
# load(paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments_pulled.RData"))

found_in_panama = function(assessment) {
  locations = assessment$locations %>%
    filter(code == "PA")
}

loc_panama = map_dfr(a_data, function(x) {
  x$location %>%
    filter(code == "PA") %>%
    mutate(assessment_id = x$assessment_id$value,
           seasonality = as.character(seasonality),
           formerlyBred = as.character(formerlyBred))
}) %>%
  mutate(iucn_panama = TRUE)

abn_panama

gbif_iucn_panama = gbif_backbone_key %>%
  left_join(abn, by = c("genus_q", "species_q")) %>%
  left_join(loc_panama, by = "assessment_id") %>%
  mutate(iucn_panama = if_else(is.na(iucn_panama) & !is.na(assessment_id), FALSE, iucn_panama)) %>%
  select(any_of(colnames(gbif_backbone_key)),
         assessment_id,
         iucn_panama) %>%
  mutate(found_in_panama = (panama_occ_count > 0) & !(iucn_panama %in% FALSE))


hit_panama = hit_table_clean %>%
  left_join(gbif_iucn_panama %>%
              select(genus_q,
                     species_q,
                     panama_occ_count,
                     iucn_panama,
                     found_in_panama), by = c("genus" = "genus_q", "species_simple" = "species_q")) %>%
  filter(found_in_panama | is.na(found_in_panama)) %>%
  select(any_of(colnames(hit_table_clean)),
         panama_occ_count,
         found_in_panama)

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

consensus_table = hit_table_clean %>%
  filter(FALSE) %>%
  select(asv_id,
         asv_total_count,
         sample_count,
         samples,
         sequence,
         method,
         class,
         order,
         family,
         genus,
         species,
         species_simple,
         scientific_name)

find_consensus = function(hit_table, taxonomic_level, identity_th = 0, overlap_th = 0, consensus_ratio = 0.5, min_consensus_hits = 1, tie_handling = "error") {
  
  if (!(tie_handling %in% c("error", "pass", "drop"))) {
    stop('Invalid value for tie_handling. Must be in c("error", "pass", "drop")')
  }
  
  hit_table_filtered = hit_table %>%
    filter(p_identical >= identity_th,
           p_overlap >= overlap_th)
  
  if (consensus_ratio < 0.5) {
    warning("Using a consensus < 0.5 may result in multiple consensus classification. In such cases, the consensus will be selected with the greatest consensus_ratio, p_identity, and p_overlap (in that order). In the case of a true tie, only one is returned.")
  }
  
  if (!(taxonomic_level %in% c("class", "order", "family", "genus", "species"))) {
    stop('taxonnomic_level not recognized, must be one of c("class", "order", "family", "genus", "species")')
  }
  
  taxon_groups_list = list(
    species = c("asv_id", "method", "class", "order", "family", "genus", "species"),
    genus = c("asv_id", "method", "class", "order", "family", "genus"),
    family = c("asv_id", "method", "class", "order", "family"),
    order = c("asv_id", "method", "class", "order"),
    class = c("asv_id", "method", "class")
  )
  
  taxon_group = taxon_groups_list[[taxonomic_level]]
  
  consensus_table = hit_table_filtered %>%
    group_by(asv_id, method) %>%
    mutate(hits_n = n()) %>%
    ungroup() %>%
    group_by_at(taxon_group) %>%
    mutate(taxon_hits_n = n()) %>%
    ungroup() %>%
    mutate(consensus_value = taxon_hits_n / hits_n) %>%
    filter(taxon_hits_n >= min_consensus_hits,
           consensus_value >= consensus_ratio) %>%
    arrange(desc(asv_total_count),
            asv_id,
            method,
            desc(consensus_value),
            desc(p_identical),
            (p_overlap)) %>%
    select(all_of(c(taxon_group,
             "consensus_value",
             "p_identical",
             "p_overlap"))) %>%
    distinct() %>%
    group_by_at(taxon_group) %>%
    slice_max(order_by = consensus_value, n = 1, with_ties = TRUE) %>%
    slice_max(order_by = p_identical, n = 1, with_ties = TRUE) %>%
    slice_max(order_by = p_overlap, n = 1, with_ties = TRUE) %>%
    mutate(consensus_n = n()) %>%
    ungroup()
  
  tie_ids = consensus_table %>%
    filter(consensus_n > 1) %>%
    select(asv_id) %>%
    distinct() %>%
    pull(asv_id)
  
  if (length(tie_ids) > 0) {
    if (tie_handling == "error") {
      stop(paste0("Ties found for the following asv_ids:\n\t", paste(collapse = "\n\t"), "\nAborting."))
    } else if (tie_handling == "pass") {
      warning(paste0("Ties found for the following asv_ids:\n\t", paste(collapse = "\n\t"), "\nPassing to output."))
    } else if (tie_handling == "drop") {
      warning(paste0("Ties found for the following asv_ids:\n\t", paste(collapse = "\n\t"), "\nDropping from output."))
      consensus_table = consensus_table %>%
        filter(consensus <= 1)
    }
  }
  
  return(consensus_table)
}

my_first_consensus = find_consensus(hit_panama, "species")

# here we need a function which takes in percentIdentical, min hits, and taxonomic level and generates a consensus, to run an itterative consensus across various criteria.

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
