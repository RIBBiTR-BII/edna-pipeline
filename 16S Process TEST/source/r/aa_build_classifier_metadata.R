library(yaml)
library(lubridate)
library(here)

# define emtadata
metadata = list()
metadata$classifier_name = "NCBI-Vertebrata16S"
metadata$classifier_script = "./source/shell/fetch-ncbi-seq-classifier-vertebrata16s.sh"
metadata$claddifier_build_timestamp = format(now(), '%Y-%m-%d %H:%M:%S')
# write to fiel
write_yaml(metadata, here("vertebrata16s-classification", "classifier_metadata.yml"))