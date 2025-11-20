# Title: Process taxonomic hits to identify species likely found in Panama
# Created by: Cob Staines (cobstainesconsulting@gmail.com)

## Load in Libraries (install as needed)
library(tidyverse)
library(iucnredlist)
library(rgbif)
library(yaml)

# config
setwd("16S Process TEST")
env_config_path = "runs/2025-11-07_panama/output/metadata/config.yml"
iucn_api_token = Sys.getenv("iucn_token")

# read in config file
config = read_yaml(env_config_path)

gbif_bool = config$geography$run_gbif_lookup
iucn_bool = config$geography$run_iucn_lookup
country_codes = config$geography$country_codes

taxa_raw = read_csv(paste0(config$run$runDir, "/output/", config$run$name, "_taxonomy_table_filtered.csv"))

taxa_clean = taxa_raw %>%
  mutate(
    species_simple = gsub(" .*", "", species),
    species_simple = if_else(species_simple %in% c("sp.",
                                                   "aff.",
                                                   "cf."),
                             NA,
                             species_simple)
  )

# copy for output
taxa_out = taxa_clean

# get unique for querying
taxa_unique = taxa_clean %>%
  select(class, order, family, genus, species_simple) %>%
  distinct() %>%
  arrange(class, order, family, genus, species_simple)

# pull from gbif
if (gbif_bool) {
  cat("Query GBIF for occurences in target countries (families, genuses & species)...\n")
  
  ## species
  cat("\tLooking up species IDs\n")
  
  gbif_query_species = taxa_unique %>%
    mutate(scientificName = case_when(
      !is.na(genus) & !is.na(species_simple) ~ paste0(genus, " ", species_simple),
      !is.na(genus) & is.na(species_simple) ~ genus,
      .default = NA)) %>%
    rename(species = species_simple) %>%
    filter(!is.na(species)) %>%
    distinct()
  
  gbif_backbone_species = name_backbone_checklist(gbif_query_species)
  
  rownames(gbif_query_species) = NULL
  rownames(gbif_backbone_species) = NULL
  
  gbif_backbone_key_species = bind_cols(gbif_query_species %>%
                                          rename(class_q = class,
                                                 order_q = order,
                                                 family_q = family,
                                                 genus_q = genus,
                                                 species_q = species,
                                                 scientific_name_q = scientificName),
                                        gbif_backbone_species) %>%
    mutate(local_occ_count = NA) %>%
    filter(rank == "SPECIES")
  
  ## genus
  cat("\tLooking up genus IDs\n")
  
  gbif_query_genus = taxa_unique %>%
    select(-species_simple) %>%
    filter(!is.na(genus)) %>%
    distinct() %>%
    mutate(scientificName = genus)
  
  gbif_backbone_genus = name_backbone_checklist(gbif_query_genus)
  
  rownames(gbif_query_genus) = NULL
  rownames(gbif_backbone_genus) = NULL
  
  gbif_backbone_key_genus = bind_cols(gbif_query_genus %>%
                                        rename(class_q = class,
                                               order_q = order,
                                               family_q = family,
                                               genus_q = genus,
                                               scientific_name_q = scientificName),
                                      gbif_backbone_genus) %>%
    mutate(local_occ_count = NA) %>%
    filter(rank == "GENUS")
  
  ## family
  cat("\tLooking up family IDs\n")
  
  gbif_query_family = taxa_unique %>%
    select(-genus,
           -species_simple) %>%
    filter(!is.na(family)) %>%
    distinct() %>%
    mutate(scientificName = family)
  
  gbif_backbone_family = name_backbone_checklist(gbif_query_family)
  
  rownames(gbif_query_family) = NULL
  rownames(gbif_backbone_family) = NULL
  
  gbif_backbone_key_family = bind_cols(gbif_query_family %>%
                                         rename(class_q = class,
                                                order_q = order,
                                                family_q = family,
                                                scientific_name_q = scientificName),
                                       gbif_backbone_family) %>%
    mutate(local_occ_count = NA) %>%
    filter(rank == "FAMILY")
  
  # row bind
  gbif_backbone_key = bind_rows(gbif_backbone_key_species,
                                gbif_backbone_key_genus,
                                gbif_backbone_key_family)
  
  occ_search
  
  cat("\tQuerying local taxa occurences\n")
  pb = txtProgressBar(min = 1, max = nrow(gbif_backbone_key), style = 3)
  for (ii in 1:nrow(gbif_backbone_key)) {
    if (!is.na(gbif_backbone_key$usageKey[ii])) {
      gbif_backbone_key$local_occ_count[ii] = occ_count(taxonKey = gbif_backbone_key$usageKey[ii],
                                                         country = country_codes)
      Sys.sleep(0.1)
      setTxtProgressBar(pb, ii)
    }
  }
  close(pb)
  cat("Done.\n")
  
  # join with taxa_out
  taxa_out = taxa_out %>%
    left_join(gbif_backbone_key %>%
                filter(rank == "FAMILY") %>%
                rename(local_gbif_count_family = local_occ_count,
                       gbif_familyKey = familyKey) %>%
                select(class_q,
                       order_q,
                       family_q,
                       gbif_familyKey,
                       local_gbif_count_family),
              by = c("class" = "class_q",
                     "order" = "order_q",
                     "family" = "family_q")) %>%
    left_join(gbif_backbone_key %>%
                filter(rank == "GENUS") %>%
                rename(local_gbif_count_genus = local_occ_count,
                       gbif_genusKey = genusKey) %>%
                select(class_q,
                       order_q,
                       family_q,
                       genus_q,
                       gbif_genusKey,
                       local_gbif_count_genus),
              by = c("class" = "class_q",
                     "order" = "order_q",
                     "family" = "family_q",
                     "genus" = "genus_q")) %>%
    left_join(gbif_backbone_key %>%
                filter(rank == "SPECIES") %>%
                rename(local_gbif_count_species = local_occ_count,
                       gbif_speciesKey = speciesKey) %>%
                select(class_q,
                       order_q,
                       family_q,
                       genus_q,
                       species_q,
                       gbif_speciesKey,
                       local_gbif_count_species),
              by = c("class" = "class_q",
                     "order" = "order_q",
                     "family" = "family_q",
                     "genus" = "genus_q",
                     "species_simple" = "species_q"))
  
  
  # write_csv(gbif_backbone_key, paste0(config$run$runDir, "/output/", config$run$name, "_GBIF_panama_occurence.csv"))
  # gbif_backbone_key = read_csv(paste0(config$run$runDir, "/output/", config$run$name, "_GBIF_panama_occurence.csv"))
}

if (iucn_bool) {
  cat("Query IUCN for occurences in target countries (species only)...\n")
  
  # cross-reference with IUCN
  api = init_api(iucn_api_token)
  
  iucn_query_species = taxa_unique %>%
    rename(species = species_simple) %>%
    filter(!is.na(genus) & !is.na(species)) %>%
    distinct()
  
  # init empty assessment df
  iucn_assessment = assessments_by_name(api, genus = iucn_query_species$genus[1], species = iucn_query_species$species[1]) %>%
    mutate(genus_q = iucn_query_species$genus[1],
           species_q = iucn_query_species$species[1]) %>%
    filter(FALSE) %>%
    select(genus_q,
           species_q,
           everything())
  
  cat("\tLooking up species IDs\n")
  pb = txtProgressBar(min = 1, max = nrow(iucn_query_species), style = 3)
  for (ii in 1:nrow(iucn_query_species)) {
    ii_genus = iucn_query_species$genus[ii]
    ii_species = iucn_query_species$species[ii]
    setTxtProgressBar(pb, ii)
    
    # pull data
    new_row <- tryCatch({
      assessments_by_name(api, genus = ii_genus, species = ii_species) %>%
        mutate(genus_q = ii_genus,
               species_q = ii_species) %>%
        slice_max(order_by = year_published, n = 1)
    }, error = function(e) {
      # fail silently
      NULL
    })
    
    # bind
    if (!is.null(new_row) && nrow(new_row) > 0) {
      iucn_assessment = bind_rows(iucn_assessment, new_row)
    }
    # sleep to avoid rate limiting
    Sys.sleep(0.5)
  }
  close(pb)
  
 #  write_csv(iucn_assessment, paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments.csv"))
  # iucn_assessment = read_csv(paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments.csv"))
  
  cat("\tQuerying local taxa occurences\n")
  iucn_assessment_data = assessment_data_many(api,
                                              iucn_assessment$assessment_id,
                                              wait_time = 0.5)
  # saveRDS(a_data, file = paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments_pulled.RData"))
  # load(paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments_pulled.RData"))
  
  f_found_locally = function(assessment) {
    locations = assessment$locations %>%
      filter(code %in% country_codes) %>%
      mutate(assessment_id = assessment$assessment_id$value,
             seasonality = as.character(seasonality),
             formerlyBred = as.character(formerlyBred))
  }
  
  local_assessments = map_dfr(iucn_assessment_data, found_locally)
  
  iucn_assessment_local = iucn_assessment %>%
    left_join(local_assessments, by = "assessment_id") %>%
    rename(iucn_taxon_id = sis_taxon_id) %>%
    group_by(genus_q, species_q, iucn_taxon_id) %>%
    summarize(iucn_reported_locally = if_else(any(!is.na(code)), TRUE, FALSE),
              .groups = "drop")
  
  # join with taxa_out
  taxa_out = taxa_out %>%
    left_join(iucn_assessment_local, by = c("genus" = "genus_q",
                                           "species_simple" = "species_q"))
}

write.csv(taxa_out, paste0(config$run$runDir, "/output/", config$run$name, "_taxonomy_table_filtered_geography.csv"), row.names = FALSE)
