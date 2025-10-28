library(tidyverse)

# create manifest file in tabular form
tibble(`absolute-filepath` = c(list.files('sequences', recursive = TRUE))) %>%
  mutate(`sample-name` = str_match(`absolute-filepath`, "(Pe-[0-9A-Z-]*_S[0-9]*)_L[0-9]*_R[12]")[,2],
         `sample-id` = str_match(`absolute-filepath`, "(Pe-[0-9A-Z-]*)_S[0-9]*_L[0-9]*_R[12]")[,2],
         `absolute-filepath` = paste0("$PWD/sequences/", `absolute-filepath`),
         `direction` = if_else(str_detect(`absolute-filepath`, "_R1_"), "forward", 'reverse'),
         `sample-id` = if_else(`sample-id` == 0), `sample-name`, `sample-id`) %>%
  select(`sample-id`, `absolute-filepath`, `direction`) %>%
  mutate(`sample-id` = str_replace(`sample-id`, "-", "_")) %>%
  # write as a csv and place in metadata folder
  write_csv(., "metadata/manifest.csv")
