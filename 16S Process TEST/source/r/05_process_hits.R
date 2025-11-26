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
  select(-consensus) %>%
  mutate(p_identical = percent_identical / 100,
         p_overlap = seq_overlap / (s_end - s_start + 1))

taxa_geography = read.csv(paste0(config$run$runDir, "/output/", config$run$name, "_taxonomy_table_filtered_geography.csv"))

hits_vsearch_geography = hits_vsearch_clean %>%
  left_join(taxa_geography %>%
              select(taxon_id,
                     any_of(c("specificEpithet",
                              "intraspecificEpithet",
                              "gbif_count_family_local",
                              "gbif_count_genus_local",
                              "gbif_count_species_local",
                              "iucn_reported_species_local"))), by = "taxon_id") %>%
  mutate(local_species = case_when(
    iucn_reported_species_local %in% TRUE ~ TRUE,
    gbif_count_species_local > 0 ~ TRUE,
    iucn_reported_species_local %in% FALSE ~ FALSE,
    gbif_count_species_local == 0 ~ FALSE,
    gbif_count_genus_local == 0 ~ FALSE,
    gbif_count_family_local == 0 ~ FALSE,
    .default = NA))

find_consensus = function(hit_table, taxonomic_level, identity_th = 0, consensus_th = 0.5, min_consensus_hits = 1, tie_handling = "error") {
  
  if (!(tie_handling %in% c("error", "pass", "drop"))) {
    stop('Invalid value for tie_handling. Must be in c("error", "pass", "drop")')
  }
  
  hit_table_filtered = hit_table %>%
    filter(p_identical >= identity_th)
  
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
    mutate(asv_hits_n = n()) %>%
    ungroup() %>%
    group_by_at(taxon_group) %>%
    mutate(taxon_hits_n = n()) %>%
    ungroup() %>%
    mutate(consensus_value = taxon_hits_n / asv_hits_n) %>%
    filter(taxon_hits_n >= min_consensus_hits,
           consensus_value >= consensus_th) %>%
    arrange(desc(asv_total_count),
            asv_id,
            method,
            desc(consensus_value),
            desc(p_identical)) %>%
    group_by_at(c(taxon_group,
                  "consensus_value")) %>%
    summarise(n_hits = first(taxon_hits_n),
              p_identical_min = mean(p_identical),
              p_identical_max = mean(p_identical),
              p_identical_mean = mean(p_identical),
              seq_overlap_max = max(seq_overlap),
              .groups = "drop") %>%
    group_by(asv_id, method) %>%
    slice_max(order_by = consensus_value, n = 1, with_ties = TRUE) %>%
    mutate(n_hits_consensus = n()) %>%
    ungroup() %>%
    mutate(identity_th = identity_th,
           consensus_th = consensus_th,
           min_consensus_hits = min_consensus_hits,
           taxonomic_level = taxonomic_level)
  
  tie_ids = consensus_table %>%
    filter(n_hits_consensus > 1) %>%
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

hits_group_calc = function(hits_df, group_cols) {
  hits_df %>%
    group_by_at(group_cols) %>%
    summarize(
      n_hits = n(),
      p_identical_mean = mean(percent_identical) / 100,
      p_identical_max = max(percent_identical) / 100,
      p_identical_min = min(percent_identical) / 100,
      seq_overlap_max = max(seq_overlap),
      .groups = "drop"
    ) %>%
    arrange(desc(asv_total_count), asv_id, method, desc(p_identical_max)) %>%
    group_by(asv_id, method) %>%
    mutate(next_pimax = lead(p_identical_max),
           drop_to_next_pimax = p_identical_max - next_pimax,
           pimax_rank = row_number()) %>%
    ungroup() %>%
    select(all_of(group_cols),
           n_hits,
           p_identical_max,
           p_identical_mean,
           p_identical_min,
           drop_to_next_pimax,
           everything())
}

hits_vsearch_species_local = hits_group_calc(hits_vsearch_geography %>%
                                               filter(!(local_species %in% FALSE),
                                                      !is.na(specificEpithet)),
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
                                               "specificEpithet",
                                               "gbif_count_family_local",
                                               "gbif_count_genus_local",
                                               "gbif_count_species_local",
                                               "iucn_reported_species_local",
                                               "local_species"))

hits_vsearch_genus_local = hits_group_calc(hits_vsearch_geography %>%
                                             filter(gbif_count_genus_local > 20,
                                                    !is.na(genus)),
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
                                             "gbif_count_family_local",
                                             "gbif_count_genus_local"))

hits_vsearch_family_local = hits_group_calc(hits_vsearch_geography %>%
                                              filter(gbif_count_family_local > 50,
                                                     !is.na(family)),
                                            c("asv_id",
                                              "asv_total_count",
                                              "sample_count",
                                              "samples",
                                              "sequence",
                                              "method",
                                              "class",
                                              "order",
                                              "family",
                                              "gbif_count_family_local"))

hits_vsearch_order = hits_group_calc(hits_vsearch_geography %>%
                                       filter(!is.na(order)),
                                     c("asv_id",
                                       "asv_total_count",
                                       "sample_count",
                                       "samples",
                                       "sequence",
                                       "method",
                                       "class",
                                       "order"))

accept_species_100_consensus = find_consensus(hits_vsearch_geography %>%
                                                filter(!is.na(order)), "species", identity_th = 1, consensus_th = 0.5001, tie_handling = "drop") %>%
  mutate(method = "1_species_100p_consensus") %>%
  left_join(hits_vsearch_geography %>%
              select(asv_id,
                     asv_total_count,
                     sample_count,
                     samples,
                     sequence) %>%
              distinct(), by = "asv_id")

accept_species_100_unique = hits_vsearch_geography %>%
  filter(p_identical == 1,
         !is.na(species)) %>%
  group_by(asv_id,
           asv_total_count,
           sample_count,
           samples,
           sequence,
           class,
           order,
           family,
           genus,
           specificEpithet) %>%
  summarise(n_hits = n(),
            p_identical_mean = mean(percent_identical) / 100,
            p_identical_max = max(percent_identical) / 100,
            p_identical_min = min(percent_identical) / 100,
            seq_overlap_max = max(seq_overlap),
            gbif_count_family_local = first(gbif_count_family_local),
            gbif_count_genus_local = first(gbif_count_genus_local),
            gbif_count_species_local = first(gbif_count_species_local),
            iucn_reported_species_local = any(iucn_reported_species_local, na.rm = TRUE),
            local_species = any(local_species, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(method = "1_species_100p_unique") %>%
  group_by(asv_id) %>%
  mutate(tie = n() > 1) %>%
  ungroup()

accept_species_local_unambiguous = hits_vsearch_species_local %>%
  # anti_join(accept_species_100_consensus, by = "asv_id") %>%
  filter(pimax_rank == 1,
         p_identical_max >= 0.95,
         is.na(drop_to_next_pimax) | drop_to_next_pimax >= 0.03) %>%
  mutate(method = "2_species_local_95p_pdrop_03",
         species = specificEpithet)

accept_genus_local_unambiguous = hits_vsearch_genus_local %>%
  # anti_join(accept_species_100_consensus, by = "asv_id") %>%
  # anti_join(accept_species_local_unambiguous, by = "asv_id") %>%
  filter(pimax_rank == 1,
         p_identical_max >= 0.90,
         is.na(drop_to_next_pimax) | drop_to_next_pimax >= 0.03) %>%
  mutate(method = "3_genus_local_90p_pdrop_03")

accept_family_local_unambiguous = hits_vsearch_family_local %>%
  # anti_join(accept_species_100_consensus, by = "asv_id") %>%
  # anti_join(accept_species_local_unambiguous, by = "asv_id") %>%
  # anti_join(accept_genus_local_unambiguous, by = "asv_id") %>%
  filter(pimax_rank == 1,
         p_identical_max >= 0.80,
         is.na(drop_to_next_pimax) | drop_to_next_pimax >= 0.05) %>%
  mutate(method = "4_family_local_80p_pdrop_05")

accept_order_unambiguous = hits_vsearch_order %>%
  # anti_join(accept_species_100_consensus, by = "asv_id") %>%
  # anti_join(accept_species_local_unambiguous, by = "asv_id") %>%
  # anti_join(accept_genus_local_unambiguous, by = "asv_id") %>%
  # anti_join(accept_family_local_unambiguous, by = "asv_id") %>%
  filter(pimax_rank == 1,
         p_identical_max >= 0.70,
         is.na(drop_to_next_pimax) | drop_to_next_pimax >= 0.05) %>%
  mutate(method = "5_order_70p_pdrop_05")

accept_order_70_consensus = find_consensus(hits_vsearch_geography,
                                           # anti_join(accept_species_100_consensus, by = "asv_id") %>%
                                           # anti_join(accept_species_local_unambiguous, by = "asv_id") %>%
                                           # anti_join(accept_genus_local_unambiguous, by = "asv_id") %>%
                                           # anti_join(accept_family_local_unambiguous, by = "asv_id") %>%
                                           # anti_join(accept_order_unambiguous, by = "asv_id"),
                                           "order", identity_th = 0.7, consensus_th = 0.5001, tie_handling = "drop") %>%
  mutate(method = "6_order_70p_consensus") %>%
  left_join(hits_vsearch_geography %>%
              select(asv_id,
                     asv_total_count,
                     sample_count,
                     samples,
                     sequence) %>%
              distinct(), by = "asv_id")

accept_all = bind_rows(accept_species_local_unambiguous,
                       accept_genus_local_unambiguous,
                       accept_family_local_unambiguous,
                       accept_order_unambiguous,
                       # accept_species_100_consensus,
                       accept_species_100_unique,
                       accept_order_70_consensus) %>%
  mutate(tie = if_else(is.na(tie), FALSE, tie)) %>%
  select(-drop_to_next_pimax,
         -next_pimax,
         -pimax_rank) %>%
  arrange(desc(asv_total_count), asv_id, method, desc(n_hits), desc(p_identical_max))

asv_count = length(unique(accept_all$asv_id))

remainder = hits_vsearch_geography %>%
  anti_join(accept_all, by = "asv_id")

hits_accept = hits_vsearch_geography %>%
  bind_rows(accept_all) %>%
  arrange(desc(asv_total_count), asv_id, method, desc(n_hits), desc(p_identical))

write.csv(hits_accept, paste0(config$run$runDir, "/output/", config$run$name, "_vsearch_hybrid_classification.csv"), row.names = FALSE)

# checks
iucn_no_gbif_yes = taxa_geography %>%
  select(-taxon_id) %>%
  distinct() %>%
  filter(iucn_reported_species_local %in% FALSE,
         gbif_count_species_local > 0) %>%
  arrange(desc(gbif_count_species_local))

iucn_yes_gbif_no = taxa_geography %>%
  select(-taxon_id) %>%
  distinct() %>%
  filter(iucn_reported_species_local %in% TRUE,
         gbif_count_species_local == 0) %>%
  arrange(desc(gbif_count_species_local))

iucn_no_gbif_no = taxa_geography %>%
  select(-taxon_id) %>%
  distinct() %>%
  filter(iucn_reported_species_local %in% FALSE,
         gbif_count_species_local == 0) %>%
  arrange(desc(gbif_count_species_local))

nonlocal_100p = accept_species_100_consensus %>%
  left_join(taxa_geography %>%
              select(-taxon_id) %>%
              distinct(), by = c("class", "order", "family", "genus", "species")) %>%
  filter(iucn_reported_species_local %in% FALSE & gbif_count_species_local == 0)

short_seq = hits_accept %>%
  filter(asv_id %in% (hits_accept %>%
                        filter(seq_overlap_max < 200) %>%
                        pull(asv_id) %>%
                        unique()))

hits_accept_100p_ties = hits_accept %>%
  filter(asv_id %in% (hits_accept %>%
                        filter(tie) %>%
                        pull(asv_id) %>%
                        unique()))

write.csv(hits_accept_100p_ties, paste0(config$run$runDir, "/output/", config$run$name, "_vsearch_hybrid_classification_ties.csv"), row.names = FALSE)

accept_out = accept_all %>%
  group_by(asv_id) %>%
  slice_min(method, n = 1, with_ties = TRUE) %>%
  ungroup() %>%
  mutate(species = if_else(is.na(species), specificEpithet, species),
         taxonomic_level = case_when(grepl("species", method) ~ "species",
                                     grepl("genus", method) ~ "genus",
                                     grepl("family", method) ~ "family",
                                     grepl("order", method) ~ "order")) %>%
  select(asv_id,
         method,
         taxonomic_level,
         class,
         order,
         family,
         genus,
         species,
         n_hits,
         p_identical_max,
         seq_overlap_max,
         consensus_value,
         tie)

write.csv(accept_out, paste0(config$run$runDir, "/output/", config$run$name, "_asv_classified.csv"), row.names = FALSE)
