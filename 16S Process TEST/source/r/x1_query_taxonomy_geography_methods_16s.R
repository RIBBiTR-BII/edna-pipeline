# Title: Process taxonomic hits to identify species likely found in Panama
# Created by: Cob Staines (cobstainesconsulting@gmail.com)

## Load in Libraries (install as needed)
library(tidyverse)
library(iucnredlist)
library(rgbif)
library(yaml)

# config
# setwd("16S Process TEST")
print(getwd())
env_config_path = "runs/methods_2026-01-12/16S/output/metadata/config.yml"
iucn_api_token = Sys.getenv("iucn_token")

# read in config file
config = read_yaml(env_config_path)

gbif_bool = config$geography$run_gbif_lookup
iucn_bool = config$geography$run_iucn_lookup
country_codes = config$geography$country_codes


taxa_raw = read_csv(paste0(config$run$runDir, "/output/", config$run$name, "_taxonomy_table_filtered.csv"))

taxa_clean = taxa_raw %>%
  mutate(
    scientificName = case_when(
      !is.na(genus) & !is.na(species) ~ paste0(genus, " ", species),
      !is.na(genus) & is.na(species) ~ genus,
      .default = NA),
    species_simple = species,
    species_simple = if_else(grepl("^sp\\.", species_simple), NA, species_simple),
    species_simple = gsub(" x .*", "", species_simple),
    species_simple = gsub("^cf\\. *", "", species_simple),
    species_simple = gsub("^aff\\. *", "", species_simple),
    specificEpithet = gsub(" .*", "", species_simple),
    intraspecificEpithet = gsub("^[^ ]+ ?", "", species_simple),
    intraspecificEpithet = if_else(intraspecificEpithet == "", NA, intraspecificEpithet))

# copy for output
taxa_out = taxa_clean %>%
  mutate(country_codes = paste(country_codes, collapse = ", "))

# get unique for querying
taxa_unique = taxa_clean %>%
  select(-taxon_id,
         -kingdom,
         -phylum) %>%
  distinct() %>%
  arrange(class, order, family, genus, species)

# pull from gbif
if (gbif_bool) {
  cat("\nQuery GBIF for occurences in target countries (families, genuses & species)...\n")
  
  ## species
  cat("\tLooking up species IDs\n")
  
  # by full scientific name
  gbif_query_species_sn = taxa_unique %>%
    filter(!is.na(species)) %>%
    select(class,
           order,
           family,
           scientificName) %>%
    distinct()
  
  gbif_backbone_species_sn = name_backbone_checklist(gbif_query_species_sn, bucket_size = 100)
  
  rownames(gbif_query_species_sn) = NULL
  rownames(gbif_backbone_species_sn) = NULL
  
  gbif_backbone_key_species_sn = bind_cols(gbif_query_species_sn %>%
                                             rename(class_q = class,
                                                    order_q = order,
                                                    family_q = family,
                                                    scientific_name_q = scientificName),
                                           gbif_backbone_species_sn) %>%
    mutate(local_occ_count = NA) %>%
    filter(rank %in% c("SPECIES", "SUBSPECIES"))
  
  # by specific epithet, for those which returned NA
  gbif_query_species_e = taxa_unique %>%
    filter(!is.na(specificEpithet)) %>%
    select(class,
           order,
           family,
           genus,
           specificEpithet,
           scientificName) %>%
    anti_join(gbif_backbone_key_species_sn, by = c("class" = "class_q",
                                                   "order" = "order_q",
                                                   "family" = "family_q",
                                                   "scientificName" = "scientific_name_q")) %>%
    select(-scientificName) %>%
    distinct()
  
  gbif_backbone_species_e = name_backbone_checklist(gbif_query_species_e, bucket_size = 100)
  
  rownames(gbif_query_species_e) = NULL
  rownames(gbif_backbone_species_e) = NULL
  
  gbif_backbone_key_species_e = bind_cols(gbif_query_species_e %>%
                                            rename(class_q = class,
                                                   order_q = order,
                                                   family_q = family,
                                                   genus_q = genus,
                                                   specific_epithet_q = specificEpithet),
                                          gbif_backbone_species_e) %>%
    mutate(local_occ_count = NA) %>%
    filter(rank %in% c("SPECIES", "SUBSPECIES"))
  
  # combine specific results
  gbif_backbone_key_species = bind_rows(gbif_backbone_key_species_sn,
                                        gbif_backbone_key_species_e)
  
  ## genus
  cat("\tLooking up genus IDs\n")
  
  gbif_query_genus = taxa_unique %>%
    select(class,
           order,
           family,
           genus) %>%
    filter(!is.na(genus)) %>%
    distinct()
  
  gbif_backbone_genus = name_backbone_checklist(gbif_query_genus, bucket_size = 100)
  
  rownames(gbif_query_genus) = NULL
  rownames(gbif_backbone_genus) = NULL
  
  gbif_backbone_key_genus = bind_cols(gbif_query_genus %>%
                                        rename(class_q = class,
                                               order_q = order,
                                               family_q = family,
                                               genus_q = genus),
                                      gbif_backbone_genus) %>%
    mutate(local_occ_count = NA) %>%
    filter(rank == "GENUS")
  
  ## family
  cat("\tLooking up family IDs\n")
  
  gbif_query_family = taxa_unique %>%
    select(class,
           order,
           family) %>%
    filter(!is.na(family)) %>%
    distinct()
  
  gbif_backbone_family = name_backbone_checklist(gbif_query_family, bucket_size = 100)
  
  rownames(gbif_query_family) = NULL
  rownames(gbif_backbone_family) = NULL
  
  gbif_backbone_key_family = bind_cols(gbif_query_family %>%
                                         rename(class_q = class,
                                                order_q = order,
                                                family_q = family),
                                       gbif_backbone_family) %>%
    mutate(local_occ_count = NA) %>%
    filter(rank == "FAMILY")
  
  # row bind
  gbif_backbone_key = bind_rows(gbif_backbone_key_species,
                                gbif_backbone_key_genus,
                                gbif_backbone_key_family)
  
  gbif_backbone_key_pa = gbif_backbone_key
  gbif_backbone_key_us = gbif_backbone_key
  gbif_backbone_key_pe = gbif_backbone_key
  
  cat("\tQuerying local taxa occurences\n")
  
  pb = txtProgressBar(min = 1, max = nrow(gbif_backbone_key), style = 3)
  for (ii in 1:nrow(gbif_backbone_key)) {
    if (!is.na(gbif_backbone_key$usageKey[ii])) {
      gbif_backbone_key_pa$local_occ_count[ii] = occ_count(taxonKey = gbif_backbone_key$usageKey[ii],
                                                        country = "PA")
      gbif_backbone_key_us$local_occ_count[ii] = occ_count(taxonKey = gbif_backbone_key$usageKey[ii],
                                                           country = "US")
      gbif_backbone_key_pe$local_occ_count[ii] = occ_count(taxonKey = gbif_backbone_key$usageKey[ii],
                                                        country = "US",
                                                        stateProvince = "Pennsylvania")
      Sys.sleep(0.2)
      setTxtProgressBar(pb, ii)
    }
  }
  close(pb)
  
  gbif_backbone_key_pa = gbif_backbone_key_pa %>%
    rename(occ_count_pa = local_occ_count)
  gbif_backbone_key_us = gbif_backbone_key_us %>%
    rename(occ_count_us = local_occ_count)
  gbif_backbone_key_pe = gbif_backbone_key_pe %>%
    rename(occ_count_pe = local_occ_count)
  
  gbif_backbone_key = bind_cols(gbif_backbone_key_pa,
                                gbif_backbone_key_us %>%
                                  select(occ_count_us),
                                gbif_backbone_key_pe %>%
                                  select(occ_count_pe))
  
  # join with taxa_out
  taxa_out = taxa_out %>%
    left_join(gbif_backbone_key %>%
                filter(rank == "FAMILY") %>%
                rename(gbif_count_family_pa = occ_count_pa,
                       gbif_count_family_us = occ_count_us,
                       gbif_count_family_pe = occ_count_pe,
                       gbif_familyKey = familyKey) %>%
                select(class_q,
                       order_q,
                       family_q,
                       gbif_familyKey,
                       gbif_count_family_pa,
                       gbif_count_family_us,
                       gbif_count_family_pe),
              by = c("class" = "class_q",
                     "order" = "order_q",
                     "family" = "family_q")) %>%
    left_join(gbif_backbone_key %>%
                filter(rank == "GENUS") %>%
                rename(gbif_count_genus_pa = occ_count_pa,
                       gbif_count_genus_us = occ_count_us,
                       gbif_count_genus_pe = occ_count_pe,
                       gbif_genusKey = genusKey) %>%
                select(class_q,
                       order_q,
                       family_q,
                       genus_q,
                       gbif_genusKey,
                       gbif_count_genus_pa,
                       gbif_count_genus_us,
                       gbif_count_genus_pe),
              by = c("class" = "class_q",
                     "order" = "order_q",
                     "family" = "family_q",
                     "genus" = "genus_q")) %>%
    left_join(gbif_backbone_key %>%
                filter(rank %in% c("SPECIES", "SUBSPECIES"),
                       !is.na(scientific_name_q)) %>%
                rename(gbif_count_species_pa_s = occ_count_pa,
                       gbif_count_species_us_s = occ_count_us,
                       gbif_count_species_pe_s = occ_count_pe,
                       gbif_speciesKey_s = speciesKey) %>%
                select(class_q,
                       order_q,
                       family_q,
                       scientific_name_q,
                       gbif_speciesKey_s,
                       gbif_count_species_pa_s,
                       gbif_count_species_us_s,
                       gbif_count_species_pe_s),
              by = c("class" = "class_q",
                     "order" = "order_q",
                     "family" = "family_q",
                     "scientificName" = "scientific_name_q")) %>%
    left_join(gbif_backbone_key %>%
                filter(rank %in% c("SPECIES", "SUBSPECIES"),
                       !is.na(specific_epithet_q)) %>%
                rename(gbif_count_species_pa_e = occ_count_pa,
                       gbif_count_species_us_e = occ_count_us,
                       gbif_count_species_pe_e = occ_count_pe,
                       gbif_speciesKey_e = speciesKey) %>%
                select(class_q,
                       order_q,
                       family_q,
                       genus_q,
                       specific_epithet_q,
                       gbif_speciesKey_e,
                       gbif_count_species_pa_e,
                       gbif_count_species_us_e,
                       gbif_count_species_pe_e) %>%
                distinct(),
              by = c("class" = "class_q",
                     "order" = "order_q",
                     "family" = "family_q",
                     "genus" = "genus_q",
                     "specificEpithet" = "specific_epithet_q")) %>%
    mutate(gbif_speciesKey = coalesce(gbif_speciesKey_s, gbif_speciesKey_e),
           gbif_count_species_pa = coalesce(gbif_count_species_pa_s, gbif_count_species_pa_e),
           gbif_count_species_us = coalesce(gbif_count_species_us_s, gbif_count_species_us_e),
           gbif_count_species_pe = coalesce(gbif_count_species_pe_s, gbif_count_species_pe_e),
           # correct errant parents with counts less than children
           gbif_count_genus_pa = if_else((!is.na(gbif_count_species_pa) & is.na(gbif_count_genus_pa)) | 
                                              gbif_count_species_pa > gbif_count_genus_pa,
                                            gbif_count_species_pa,
                                            gbif_count_genus_pa),
           gbif_count_family_pa = if_else((!is.na(gbif_count_genus_pa) & is.na(gbif_count_family_pa)) | 
                                               gbif_count_genus_pa > gbif_count_family_pa,
                                             gbif_count_genus_pa,
                                             gbif_count_family_pa),
           gbif_count_genus_us = if_else((!is.na(gbif_count_species_us) & is.na(gbif_count_genus_us)) | 
                                           gbif_count_species_us > gbif_count_genus_us,
                                         gbif_count_species_us,
                                         gbif_count_genus_us),
           gbif_count_family_us = if_else((!is.na(gbif_count_genus_us) & is.na(gbif_count_family_us)) | 
                                            gbif_count_genus_us > gbif_count_family_us,
                                          gbif_count_genus_us,
                                          gbif_count_family_us),
           gbif_count_genus_pe = if_else((!is.na(gbif_count_species_pe) & is.na(gbif_count_genus_pe)) | 
                                           gbif_count_species_pe > gbif_count_genus_pe,
                                         gbif_count_species_pe,
                                         gbif_count_genus_pe),
           gbif_count_family_pe = if_else((!is.na(gbif_count_genus_pe) & is.na(gbif_count_family_pe)) | 
                                            gbif_count_genus_pe > gbif_count_family_pe,
                                          gbif_count_genus_pe,
                                          gbif_count_family_pe)) %>%
    select(-gbif_speciesKey_s,
           -gbif_speciesKey_e,
           -gbif_count_species_pa_s,
           -gbif_count_species_pa_e,
           -gbif_count_species_us_s,
           -gbif_count_species_us_e,
           -gbif_count_species_pe_s,
           -gbif_count_species_pe_e)
  
  gbif_count_family_pacat("Done.\n")
}

if (iucn_bool) {
  cat("\nQuery IUCN for occurences in target countries (species only)...\n")
  
  api = init_api(iucn_api_token)
  
  # define assessment function
  fetch_iucn_assessments = function(api, query_df, sleep_s = 0.5) {
    # init empty assessment df
    iucn_assessment = assessments_by_name(api, genus = query_df$genus[1], species = query_df$species[1]) %>%
      mutate(genus_q = query_df$genus[1],
             species_q = query_df$species[1]) %>%
      filter(FALSE) %>%
      select(genus_q,
             species_q,
             everything())
    
    cat("\tLooking up species IDs\n")
    pb = txtProgressBar(min = 1, max = nrow(query_df), style = 3)
    for (ii in 1:nrow(query_df)) {
      ii_genus = query_df$genus[ii]
      ii_species = query_df$species[ii]
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
      Sys.sleep(sleep_s)
    }
    close(pb)
    
    return(iucn_assessment)
  }
  
  
  ## species
  iucn_query_species = taxa_unique %>%
    select(genus, species) %>%
    filter(!is.na(genus) & !is.na(species)) %>%
    distinct()
  
  iucn_assessment_species = fetch_iucn_assessments(api, iucn_query_species)
  
  ## specific epithet for species returned as null
  iucn_query_species_e = taxa_unique %>%
    anti_join(iucn_assessment_species,
              by = c("genus" = "genus_q",
                     "species" = "species_q")) %>%
    filter(!is.na(genus) & !is.na(specificEpithet),
           species != specificEpithet) %>%
    select(genus, specificEpithet) %>%
    distinct() %>%
    rename(species = specificEpithet)
  
  iucn_assessment_species_e = fetch_iucn_assessments(api, iucn_query_species_e) %>%
    rename(specific_epithet_q = species_q)
  
  # combine specific results
  iucn_assessment_raw = bind_rows(iucn_assessment_species,
                                  iucn_assessment_species_e) %>%
    select(genus_q,
           species_q,
           specific_epithet_q,
           everything())
  
  # get rid of duplicates
  iucn_assessment = iucn_assessment_raw %>%
    distinct() %>%
    mutate(scopes_priority = case_match(scopes_description_en,
                                        "Global" ~ "1_Global",
                                        "Gulf of Mexico" ~ "2_Gulf of Mexico",
                                        .default = paste0("3_", scopes_description_en))) %>%
    group_by(genus_q, species_q, specific_epithet_q) %>%
    arrange(desc(latest), scopes_priority, by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(-scopes_priority)
  
  # fetch assessment data
  cat("\tQuerying local taxa occurences\n")
  iucn_assessment_data = assessment_data_many(api,
                                              iucn_assessment$assessment_id,
                                              wait_time = 0.5)
  # saveRDS(a_data, file = paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments_pulled.RData"))
  # load(paste0(config$run$runDir, "/output/", config$run$name, "_IUCN_species_assessments_pulled.RData"))
  
  f_found_locally = function(assessment) {
    locations = assessment$locations %>%
      filter(code %in% c("PA", "US", "PEN-OO")) %>%
      mutate(assessment_id = assessment$assessment_id$value,
             seasonality = as.character(seasonality),
             formerlyBred = as.character(formerlyBred))
  }
  
  local_assessments = map_dfr(iucn_assessment_data, f_found_locally)
  
  local_assessments_unique = local_assessments %>%
    select(assessment_id, code) %>%
    mutate(present = TRUE) %>%
    pivot_wider(names_from = "code",
                values_from = "present") %>%
    rename(iucn_assessment_pa = PA,
           iucn_assessment_us = US,
           iucn_assessment_pe = `PEN-OO`,
           ) %>%
    mutate(iucn_assessment_pa = if_else(is.na(iucn_assessment_pa), FALSE, iucn_assessment_pa),
           iucn_assessment_us = if_else(is.na(iucn_assessment_us), FALSE, iucn_assessment_us),
           iucn_assessment_pe = if_else(is.na(iucn_assessment_pe), FALSE, iucn_assessment_pe))
  
  iucn_assessment_local = iucn_assessment %>%
    left_join(local_assessments_unique, by = "assessment_id") %>%
    rename(iucn_taxon_id = sis_taxon_id)
  
  iucn_assessment_local_species_s = iucn_assessment_local %>%
    filter(!is.na(species_q)) %>%
    select(-specific_epithet_q) %>%
    rename(iucn_taxon_id_s = iucn_taxon_id,
           iucn_assessment_pa_s = iucn_assessment_pa,
           iucn_assessment_us_s = iucn_assessment_us,
           iucn_assessment_pe_s = iucn_assessment_pe) %>%
    select(genus_q,
           species_q,
           iucn_taxon_id_s,
           iucn_assessment_pa_s,
           iucn_assessment_us_s,
           iucn_assessment_pe_s)
  
  iucn_assessment_local_species_e = iucn_assessment_local %>%
    filter(!is.na(specific_epithet_q)) %>%
    select(-species_q) %>%
    rename(iucn_taxon_id_e = iucn_taxon_id,
           iucn_assessment_pa_e = iucn_assessment_pa,
           iucn_assessment_us_e = iucn_assessment_us,
           iucn_assessment_pe_e = iucn_assessment_pe) %>%
    select(genus_q,
           specific_epithet_q,
           iucn_taxon_id_e,
           iucn_assessment_pa_e,
           iucn_assessment_us_e,
           iucn_assessment_pe_e)
  
  # join with taxa_out
  taxa_out = taxa_out %>%
    left_join(iucn_assessment_local_species_s, by = c("genus" = "genus_q", "species_simple" = "species_q")) %>%
    left_join(iucn_assessment_local_species_e, by = c("genus" = "genus_q", "specificEpithet" = "specific_epithet_q")) %>%
    mutate(iucn_taxon_id = coalesce(iucn_taxon_id_s, iucn_taxon_id_e),
           iucn_assessment_pa = coalesce(iucn_assessment_pa_s, iucn_assessment_pa_e),
           iucn_assessment_us = coalesce(iucn_assessment_us_s, iucn_assessment_us_e),
           iucn_assessment_pe = coalesce(iucn_assessment_pe_s, iucn_assessment_pe_e)) %>%
    select(-iucn_taxon_id_s,
           -iucn_taxon_id_e,
           -iucn_assessment_pa_s,
           -iucn_assessment_pa_e,
           -iucn_assessment_us_s,
           -iucn_assessment_us_e,
           -iucn_assessment_pe_s,
           -iucn_assessment_pe_e)
  
  cat("Done.\n")
}

write.csv(taxa_out, paste0(config$run$runDir, "/output/", config$run$name, "_taxonomy_table_filtered_geography_16s.csv"), row.names = FALSE)
