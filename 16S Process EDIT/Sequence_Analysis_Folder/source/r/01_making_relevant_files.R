# Title: Organizing data from NGS for further analysis. 
# Created by: Brandon Hoenig (brandonhoenig@gmail.com)
# Created on: 12 August 2024
# Edited on: 4 November 2024

## Load in Libraries (install them in you need them.)
librarian::shelf(tidyverse, helixcn/phylotools, readxl, openxlsx, here)
# library(tidyverse)
# library(phylotools)
# library(readxl)
# library(openxlsx)

## Read in the .biom feature table. 
feature_table <-
  read_tsv(here(wddir, "analysis/denoised_16S_eDNA/sample_table/feature-table.tsv"), skip = 1) %>%
    as.data.frame() %>%
    rename("ASV_ID" = "#OTU ID")

## Read in taxonomic classifications 
taxonomy_table <-
  read.table(paste0(wddir, "/analysis/denoised_16S_eDNA/classified_taxonomy/classified_taxonomy/", 
          list.files(here(wddir, 'analysis/denoised_16S_eDNA/classified_taxonomy/classified_taxonomy')), 
          "/data/taxonomy.tsv"),
          header = T, 
          sep = '\t')

# Read in DNA sequences 
sequence_table <-
  read.fasta(file = paste0(wddir, "/analysis/denoised_16S_eDNA/representative_sequences/", 
         list.files(here(wddir, 'analysis/denoised_16S_eDNA/representative_sequences')), 
         "/data/dna-sequences.fasta"))%>%
  rename(ASV_ID = seq.name, 
         sequence = seq.text)

## Write each table for documentation.  Proceed with analysis in R or in Excel. 
list_of_datasets <- list("Feature Table" = feature_table, 
                         "Taxonomy Table" = taxonomy_table, 
                         "Sequence Table" = sequence_table)

write.xlsx(list_of_datasets, file = here(wddir, "output/collated_eDNA_data.xlsx"))
