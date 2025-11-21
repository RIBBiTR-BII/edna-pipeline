# Title: Process taxonomic hits to identify species likely found in Panama
# Created by: Cob Staines (cobstainesconsulting@gmail.com)

## Load in Libraries (install them in you need them.)
library(tidyverse)
library(yaml)


# manual runs
# setwd("16S Process TEST")
env_config_path = "runs/2025-11-07_panama/output/metadata/config.yml"

# read in config file
config = read_yaml(env_config_path)

hits_raw = read_csv(paste0(config$run$runDir, "/output/", config$run$name, "_eDNA_hits_joined.csv"))

hits_vsearch_clean = hits_raw %>%
  filter(!is.na(class),
         method == "vsearch") %>%
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

taxa_geography = read.csv(paste0(config$run$runDir, "/output/", config$run$name, "_taxonomy_table_filtered_geography.csv"))

hits_vsearch_geography = hits_vsearch_clean %>%
  left_join(taxa_geography %>%
              select(taxon_id,
                     any_of(c("local_gbif_count_family",
                              "local_gbif_count_genus",
                              "local_gbif_count_species",
                              "iucn_reported_locally_species"))), by = "taxon_id")

hits_group_calc = function(hits_df, group_cols) {
  hits_df %>%
    group_by_at(group_cols) %>%
    summarize(
      hits_n = n(),
      p_identical_mean = mean(percent_identical) / 100,
      p_identical_max = max(percent_identical) / 100,
      p_identical_min = min(percent_identical) / 100,
      .groups = "drop"
    ) %>%
    arrange(desc(asv_total_count), asv_id, method, desc(p_identical_max)) %>%
    group_by(asv_id, method) %>%
    mutate(next_pimax = lead(p_identical_max),
           drop_to_next_pimax = p_identical_max - next_pimax,
           pimax_rank = row_number()) %>%
    ungroup() %>%
    select(all_of(group_cols),
           hits_n,
           p_identical_max,
           p_identical_mean,
           p_identical_min,
           drop_to_next_pimax,
           everything())
}

hits_vsearch_local = hits_vsearch_geography %>%
  mutate(local_species = case_when(
    iucn_reported_locally_species %in% TRUE ~ TRUE,
    iucn_reported_locally_species %in% FALSE ~ FALSE,
    local_gbif_count_species > 0 ~ TRUE,
    local_gbif_count_species == 0 ~ FALSE,
    local_gbif_count_genus == 0 ~ FALSE,
    local_gbif_count_family == 0 ~ FALSE,
    .default = NA)) %>%
  filter(!(local_species %in% FALSE))

hits_vsearch_local_species = hits_group_calc(hits_vsearch_local,
                                             c("asv_id",
                                               "asv_total_count",
                                               "sample_count",
                                               "samples",
                                               "sequence",
                                               "method",
                                               "class",
                                               "order",
                                               "family",
                                               "genus",
                                               "species",
                                               "species_simple",
                                               "scientific_name",
                                               "local_gbif_count_family",
                                               "local_gbif_count_genus",
                                               "local_gbif_count_species",
                                               "iucn_reported_locally_species"))

#### resume

hits_vsearch_species = hits_vsearch_geography %>%
  group_by(asv_id,
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
           scientific_name,
           local_gbif_count_family,
           local_gbif_count_genus,
           local_gbif_count_species,
           iucn_reported_locally_species) %>%
  summarize(
    hits_n = n(),
    p_identical_mean = mean(percent_identical) / 100,
    p_identical_max = max(percent_identical) / 100,
    p_identical_min = min(percent_identical) / 100,
    .groups = "drop"
  ) %>%
  arrange(desc(asv_total_count), asv_id, method, desc(p_identical_max)) %>%
  group_by(asv_id, method) %>%
  mutate(next_pimax = lead(p_identical_max),
         drop_to_next_pimax = p_identical_max - next_pimax,
         pimax_rank = row_number()) %>%
  ungroup() %>%
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
         scientific_name,
         hits_n,
         p_identical_max,
         p_identical_mean,
         p_identical_min,
         drop_to_next_pimax,
         everything()) %>%
  mutate(local_species = case_when(
    iucn_reported_locally_species %in% TRUE ~ TRUE,
    iucn_reported_locally_species %in% FALSE ~ FALSE,
    local_gbif_count_species > 0 ~ TRUE,
    local_gbif_count_species == 0 ~ FALSE,
    local_gbif_count_genus == 0 ~ FALSE,
    local_gbif_count_family == 0 ~ FALSE,
    .default = NA))

# # unanimous
# unanimous_species = hits_vsearch_species %>%
#   group_by(asv_id,
#            asv_total_count,
#            sample_count,
#            samples,
#            sequence,
#            method) %>%
#   mutate(taxa_n = n()) %>%
#   ungroup() %>%
#   filter(taxa_n == 1)
# 
# accept_unanimous_species_local = unanimous_species %>%
#   filter(local_species %in% TRUE)
# 
# check_unanimous_species_nonlocal = unanimous_species %>%
#   filter(!(local_species %in% TRUE))

# unambiguous species
unambiguous_species = hits_vsearch_species %>%
  filter(pimax_rank == 1,
         p_identical_max >= 0.95,
         is.na(drop_to_next_pimax) | drop_to_next_pimax >= 0.03)

accept_unambiguous_species_local = unambiguous_species %>%
  filter(local_species %in% TRUE)

check_unambiguous_species_nonlocal = unambiguous_species %>%
  filter(!(local_species %in% TRUE))

remainder = hits_vsearch_species %>%
  anti_join(unambiguous_species, by = c("asv_id", "method"))

# unambiguous genus
hits_vsearch_genus = hits_vsearch_geography %>%
  anti_join(unambiguous_species, by = c("asv_id", "method")) %>%
  group_by(asv_id,
           asv_total_count,
           sample_count,
           samples,
           sequence,
           method,
           class,
           order,
           family,
           genus,
           local_gbif_count_family,
           local_gbif_count_genus) %>%
  summarize(
    hits_n = n(),
    p_identical_mean = mean(percent_identical) / 100,
    p_identical_max = max(percent_identical) / 100,
    p_identical_min = min(percent_identical) / 100,
    .groups = "drop"
  ) %>%
  arrange(desc(asv_total_count), asv_id, method, desc(p_identical_max)) %>%
  group_by(asv_id, method) %>%
  mutate(next_pimax = lead(p_identical_max),
         drop_to_next_pimax = p_identical_max - next_pimax,
         pimax_rank = row_number()) %>%
  ungroup() %>%
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
         hits_n,
         p_identical_max,
         p_identical_mean,
         p_identical_min,
         drop_to_next_pimax,
         everything()) %>%
  mutate(local_genus = local_gbif_count_genus > 50)

unambiguous_genus = hits_vsearch_genus %>%
  filter(pimax_rank == 1,
         p_identical_max >= 0.90,
         is.na(drop_to_next_pimax) | drop_to_next_pimax >= 0.03)

accept_unambiguous_genus_local = unambiguous_genus %>%
  filter(local_genus %in% TRUE)

check_unambiguous_genus_nonlocal = unambiguous_genus %>%
  filter(!(local_genus %in% TRUE))

remainder = hits_vsearch_species %>%
  anti_join(unambiguous_species, by = c("asv_id", "method")) %>%
  anti_join(unambiguous_genus, by = c("asv_id", "method"))




#############
hits_panama = hit_table_clean %>%
  left_join(gbif_iucn_panama %>%
              filter(rank == "SPECIES") %>%
              select(family_q,
                     genus_q,
                     species_q,
                     panama_occ_count,
                     iucn_panama,
                     found_in_panama) %>%
              rename(panama_occ_count_species = panama_occ_count,
                     iucn_panama_species = iucn_panama,
                     found_in_panama_species = found_in_panama), by = c("family" = "family_q", "genus" = "genus_q", "species_simple" = "species_q")) %>%
  left_join(gbif_iucn_panama %>%
              filter(rank == "GENUS") %>%
              select(family_q,
                     genus_q,
                     panama_occ_count,
                     found_in_panama) %>%
              distinct() %>%
              rename(panama_occ_count_genus = panama_occ_count,
                     found_in_panama_genus = found_in_panama), by = c("family" = "family_q", "genus" = "genus_q")) %>%
  left_join(gbif_iucn_panama %>%
              filter(rank == "FAMILY") %>%
              select(family_q,
                     panama_occ_count,
                     found_in_panama) %>%
              distinct() %>%
              rename(panama_occ_count_family = panama_occ_count,
                     found_in_panama_family = found_in_panama), by = c("family" = "family_q")) %>%
  mutate(found_in_panama_genus = if_else(found_in_panama_species, TRUE, found_in_panama_genus),
         found_in_panama_family = if_else(found_in_panama_genus, TRUE, found_in_panama_family)) %>%
  select(any_of(colnames(hit_table_clean)),
         panama_occ_count_family,
         found_in_panama_family,
         panama_occ_count_genus,
         found_in_panama_genus,
         panama_occ_count_species,
         iucn_panama_species,
         found_in_panama_species)

find_consensus = function(hit_table, taxonomic_level, identity_th = 0, overlap_th = 0, consensus_th = 0.5, min_consensus_hits = 1, tie_handling = "error") {
  
  if (!(tie_handling %in% c("error", "pass", "drop"))) {
    stop('Invalid value for tie_handling. Must be in c("error", "pass", "drop")')
  }
  
  hit_table_filtered = hit_table %>%
    filter(p_identical >= identity_th,
           p_overlap >= overlap_th)
  
  if (consensus_th <= 0.5) {
    warning("Using a consensus_th <= 0.5 may result in multiple consensus classification. Consensus classification(s) will be selected with the greatest consensus value, with ties handled according to tie_handling ('error', 'pass', 'drop')")
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
           consensus_value >= consensus_th) %>%
    arrange(desc(asv_total_count),
            asv_id,
            method,
            desc(consensus_value),
            desc(p_identical),
            (p_overlap)) %>%
    group_by_at(c(taxon_group,
                  "consensus_value")) %>%
    summarise(taxon_hits_n = first(taxon_hits_n),
              p_identical_mean = mean(p_identical),
              p_overlap_mean = mean(p_overlap),
              .groups = "drop") %>%
    group_by(asv_id, method) %>%
    slice_max(order_by = consensus_value, n = 1, with_ties = TRUE) %>%
    mutate(consensus_n = n()) %>%
    ungroup() %>%
    mutate(identity_th = identity_th,
           overlap_th = overlap_th,
           consensus_th = consensus_th,
           min_consensus_hits = min_consensus_hits,
           taxonomic_level = taxonomic_level)
  
  tie_ids = consensus_table %>%
    filter(consensus_n > 1) %>%
    select(asv_id) %>%
    distinct() %>%
    pull(asv_id)
  
  if (length(tie_ids) > 0) {
    if (tie_handling == "error") {
      stop(paste0("Ties found for the following asv_ids:\n\t", paste(tie_ids, collapse = "\n\t"), "\nAborting."))
    } else if (tie_handling == "pass") {
      warning(paste0("Ties found for the following asv_ids:\n\t", paste(tie_ids, collapse = "\n\t"), "\nPassing to output."))
    } else if (tie_handling == "drop") {
      warning(paste0("Ties found for the following asv_ids:\n\t", paste(tie_ids, collapse = "\n\t"), "\nDropping from output."))
      consensus_table = consensus_table %>%
        filter(consensus <= 1)
    }
  }
  
  return(consensus_table)
}

rolling_consensus = function(hit_table, taxonomic_level, identity_th_ub = 1, identity_th_lb = .5, idenity_th_increment = 0.05, overlap_th = 0, consensus_th = 0.5, min_consensus_hits = 1, tie_handling = "pass") {
  identity_th_seq = seq(identity_th_ub, identity_th_lb, by = -idenity_th_increment)
  
  consensus_roll = map_dfr(identity_th_seq, ~ find_consensus(hit_table, taxonomic_level, identity_th = .x, overlap_th, consensus_th, min_consensus_hits, tie_handling)) %>%
    arrange(asv_id, method, desc(identity_th)) %>%
    group_by(asv_id, method, class, order, family, genus, species, taxon_hits_n) %>%
    slice_max(identity_th, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(asv_id, method, class, order, family, genus, species) %>%
    slice_max(taxon_hits_n, n = 1) %>%
    ungroup() %>%
    group_by(asv_id, method) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    filter(n < 2) %>%
    select(-n)
}

hits_panama_grouped = hits_panama %>%
  group_by(asv_id,
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
           scientific_name,
           panama_occ_count_family,
           panama_occ_count_genus,
           panama_occ_count_species,
           iucn_panama_species) %>%
  summarize(
    hits_n = n(),
    p_identical_mean = mean(p_identical),
    p_identical_max = max(p_identical),
    .groups = "drop"
  ) %>%
  arrange(desc(asv_total_count), asv_id, method, desc(p_identical_max)) %>%
  filter(method == "vsearch") %>%
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
         scientific_name,
         hits_n,
         p_identical_max,
         p_identical_mean,
         everything())

vsearch_unanimous = hits_panama_grouped %>%
  group_by(asv_id,
           asv_total_count,
           sample_count,
           samples,
           sequence,
           method) %>%
  mutate(taxa_n = n()) %>%
  ungroup() %>%
  filter(taxa_n == 1)

consensus_table_100 = find_consensus(hits_panama, "species", identity_th = 1, overlap_th = 1, consensus_th = 0.70, tie_handling = "pass") %>%
  mutate(priority = 1)
consensus_roll_70 = rolling_consensus(hits_panama %>%
                                        filter(!(found_in_panama_species %in% FALSE)), "species", identity_th_ub = 1, identity_th_lb = .7, idenity_th_increment = 0.01, overlap_th = 0.7, consensus_th = 0.66, min_consensus_hits = 2, tie_handling = "pass") %>%
  distinct() %>%
  mutate(priority = 2)
consensus_table_70_genus = find_consensus(hits_panama %>%
                                            filter(!(found_in_panama_genus %in% FALSE)), "genus", identity_th = .70, overlap_th = .70, consensus_th = .66, tie_handling = "pass") %>%
  mutate(priority = 3)
consensus_table_60_family = find_consensus(hits_panama %>%
                                             filter(!(found_in_panama_family %in% FALSE)), "family", identity_th = .60, overlap_th = .50, consensus_th = .66, tie_handling = "pass") %>%
  mutate(priority = 4)
consensus_table_50_order = find_consensus(hits_panama, "order", identity_th = .50, overlap_th = .50, consensus_th = .5001, tie_handling = "pass") %>%
  mutate(priority = 5)


consensus_table_all = hits_panama %>%
  select(asv_id,
         asv_total_count,
         sample_count,
         samples,
         sequence) %>%
  distinct() %>%
  left_join(bind_rows(consensus_table_100,
                      consensus_roll_70,
                      consensus_table_70_genus,
                      consensus_table_60_family,
                      consensus_table_50_order), by = "asv_id") %>%
  arrange(desc(asv_total_count), asv_id, method, desc(identity_th))

consensus_table = consensus_table_all %>%
  group_by(asv_id, method) %>%
  slice_min(priority, n = 1)

consensus_table_vsearch = consensus_table %>%
  filter(method == "vsearch")

vsearch_classification_ratio = nrow(consensus_table_vsearch) / length(unique(hit_table_clean$asv_id))
blast_classification_ratio = nrow(consensus_table_blast %>%
                                    filter(!is.na(class))) / length(unique(hit_table_clean$asv_id))

consensus_table_vsearch_all = bind_rows(consensus_table_vsearch,
                                        hit_table_raw %>%
                                          select(asv_id,
                                                 asv_total_count,
                                                 sample_count,
                                                 samples,
                                                 sequence) %>%
                                          distinct() %>%
                                          mutate(method = "vsearch") %>%
                                          anti_join(consensus_table_vsearch, by = "asv_id")) %>%
  arrange(desc(asv_total_count), asv_id, method, desc(identity_th))

consensus_table_vsearch_blast = bind_rows(consensus_table_vsearch_all,
                                          peace = consensus_table_blast %>%
                                            left_join(hit_table_raw %>%
                                                        select(asv_id,
                                                               asv_total_count,
                                                               sample_count,
                                                               samples,
                                                               sequence) %>%
                                                        distinct(), by = "asv_id")) %>%
  arrange(desc(asv_total_count), asv_id, method, desc(identity_th))

no_consensus_vsearch = hit_table_raw %>%
  filter(method == "vsearch") %>%
  anti_join(consensus_table_vsearch, by = "asv_id")

write.csv(consensus_table_vsearch_all, paste0(config$run$runDir, "/output/", config$run$name, "_vsearch_consensus_classification.csv"), row.names = FALSE)
write.csv(consensus_table_vsearch_blast, paste0(config$run$runDir, "/output/", config$run$name, "_vsearch_blast_consensus_classification.csv"), row.names = FALSE)
