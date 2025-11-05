library(yaml)
library(lubridate)

# read in config file, use as basis for config
metadata = read_yaml("config.yml")
metadata$run_completion_timestamp = format(now(), '%Y-%m-%d %H:%M:%S')
write_yaml(metadata, paste0(config$`run-directory`, "/metadata/run_metadata.yml"))
