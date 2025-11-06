library(yaml)
library(lubridate)

# read in config for runDir
config = read_yaml("config.yml")
# read in run metadata
metadata = read_yaml(paste0(config$run$runDir, "/metadata/run_metadata.yml"))
# update end time
metadata$runtime$endTimestamp = format(now(), '%Y-%m-%d %H:%M:%S')
# write to file
write_yaml(metadata, paste0(config$run$runDir, "/metadata/run_metadata.yml"))
