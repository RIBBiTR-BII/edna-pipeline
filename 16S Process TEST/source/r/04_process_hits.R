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
    p_overlap = seq_overlap / q_end,
    p_identical = percent_identical / 100
  )

species_unique = hit_table_clean %>%
  select(family, genus, species_simple, scientific_name) %>%
  distinct() %>%
  arrange(family, genus, scientific_name) %>%
  mutate()

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

# pull from gbif
## species
gbif_query_species = species_unique %>%
  filter(!is.na(species_simple)) %>%
  mutate(kingdom = "Animalia")
gbif_backbone_species = name_backbone_checklist(gbif_query_species %>%
                                                  rename(scientificName = scientific_name))
rownames(gbif_query_species) = NULL
rownames(gbif_backbone_species) = NULL

gbif_backbone_key_species = bind_cols(gbif_query_species %>%
                                        select(-kingdom) %>%
                                        rename(family_q = family,
                                               genus_q = genus,
                                               species_q = species_simple,
                                               scientific_name_q = scientific_name),
                                      gbif_backbone_species) %>%
  mutate(panama_occ_count = NA)

## genus
gbif_query_genus = species_unique %>%
  select(family, genus) %>%
  filter(!is.na(genus)) %>%
  distinct() %>%
  mutate(kingdom = "Animalia")
gbif_backbone_genus = name_backbone_checklist(gbif_query_genus %>%
                                                rename(scientificName = genus))

rownames(gbif_query_genus) = NULL
rownames(gbif_backbone_genus) = NULL

gbif_backbone_key_genus = bind_cols(gbif_query_genus %>%
                                      select(-kingdom) %>%
                                      rename(family_q = family,
                                             genus_q = genus),
                                    gbif_backbone_genus) %>%
  mutate(panama_occ_count = NA)

## family
gbif_query_family = species_unique %>%
  select(family) %>%
  filter(!is.na(family)) %>%
  distinct() %>%
  mutate(kingdom = "Animalia")
gbif_backbone_family = name_backbone_checklist(gbif_query_family %>%
                                                 rename(scientificName = family))

rownames(gbif_query_family) = NULL
rownames(gbif_backbone_family) = NULL

gbif_backbone_key_family = bind_cols(gbif_query_family %>%
                                       select(-kingdom) %>%
                                       rename(family_q = family),
                                     gbif_backbone_family) %>%
  mutate(panama_occ_count = NA)

gbif_backbone_key = bind_rows(gbif_backbone_key_species,
                              gbif_backbone_key_genus,
                              gbif_backbone_key_family) %>%
  filter(rank %in% c("FAMILY", "GENUS", "SPECIES"))


print("Searching GBIF for local species occurences")
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

gbif_scarce_species = gbif_backbone_key %>%
  filter(rank == "SPECIES",
         !is.na(panama_occ_count),
         panama_occ_count < 20,)

# cross-reference with IUCN
api = init_api(Sys.getenv("iucn_token"))

# init empty df
abn = assessments_by_name(api, genus = gbif_scarce_species$genus_q[1], species = gbif_scarce_species$species_q[1]) %>%
  mutate(genus_q = gbif_scarce_species$genus_q[1],
         species_q = gbif_scarce_species$species_q[1]) %>%
  filter(FALSE) %>%
  select(genus_q,
         species_q,
         everything())

print("Cross-referencing with IUCN for GBIF-scarce species")
pb = txtProgressBar(min = 1, max = nrow(gbif_scarce_species), style = 3)
for (ii in 1:nrow(gbif_scarce_species)) {
  ii_genus = gbif_scarce_species$genus_q[ii]
  ii_species = gbif_scarce_species$species_q[ii]
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

gbif_iucn_panama = gbif_backbone_key %>%
  left_join(abn, by = c("genus_q", "species_q")) %>%
  left_join(loc_panama, by = "assessment_id") %>%
  mutate(iucn_panama = if_else(is.na(iucn_panama) & !is.na(assessment_id), FALSE, iucn_panama)) %>%
  select(any_of(colnames(gbif_backbone_key)),
         assessment_id,
         iucn_panama) %>%
  mutate(found_in_panama = case_when(
    rank == "SPECIES" ~ iucn_panama %in% TRUE | ((panama_occ_count > 0) & !(iucn_panama %in% FALSE)),
    rank == "GENUS" ~ panama_occ_count > 100,
    rank == "FAMILY" ~ panama_occ_count > 200))

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

# filter(found_in_panama | is.na(found_in_panama)) %>%
# 
# 
# hits_not_panama = hit_table_clean %>%
#   left_join(gbif_iucn_panama %>%
#               select(genus_q,
#                      species_q,
#                      panama_occ_count,
#                      iucn_panama,
#                      found_in_panama), by = c("genus" = "genus_q", "species_simple" = "species_q")) %>%
#   filter(!found_in_panama) %>%
#   select(any_of(colnames(hit_table_clean)),
#          panama_occ_count,
#          found_in_panama)

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
