library(yaml)
library(lubridate)
library(here)

# initiate classifier metadata
metadata = list()
metadata$classifier$name = "NCBI_Vertebrata16S"
metadata$classifier$buildScript = "./source/shell/fetch-ncbi-seq-classifier-vertebrata16s.sh"
metadata$classifier$buildEndTimestamp = format(now(), '%Y-%m-%d %H:%M:%S')
# write to file
write_yaml(metadata, here("vertebrata16s_classification", "classifier_metadata.yml"))