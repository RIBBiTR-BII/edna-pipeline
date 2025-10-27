library(tidyverse)
library(here)

# create manifest file in tabular form
tibble(`absolute-filepath` = c(list.files(here('16S Process TEST', '2025-10-27_sequencing', 'sequences'), recursive = TRUE))) %>%
  mutate(`sample-id` = str_replace(`absolute-filepath`, "\\_.*", "")) %>%
  separate(`sample-id`, sep = "-", into = c("seq-id", "sample-id"), extra = 'merge') %>% 
  mutate(`absolute-filepath` = paste0("$PWD/sequences/", `absolute-filepath`)) %>%
  mutate(direction = if_else(str_detect(`absolute-filepath`, "_R1_"), "forward", 'reverse')) %>%
  select(`sample-id`, `absolute-filepath`, `direction`) %>%
  mutate(`sample-id` = str_replace(`sample-id`, "-", "_")) %>%
# write as a csv and place in metadata folder
  write_csv(., "metadata/manifest.csv")
