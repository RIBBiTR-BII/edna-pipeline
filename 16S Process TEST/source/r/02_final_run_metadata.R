library(yaml)

# inherit env_config_path
args = commandArgs(trailingOnly = TRUE)
env_config_path = args[1]

# read in config for runDir
config = read_yaml(env_config_path)
# read in run metadata
metadata = read_yaml(paste0(config$run$runDir, "/metadata/run_metadata.yml"))
# update end time
end_time = Sys.time()
metadata$post$runtime$endTimestamp = format(end_time, '%Y-%m-%d %H:%M:%S %z')
# calculate duration
start_time = as.POSIXct(metadata$post$runtime$startTimestamp, format = "%Y-%m-%d %H:%M:%S %z")
metadata$post$runtime$duration = paste0(as.character(round(difftime(end_time, start_time, units = "secs"), digits = 0)), " seconds")
# update exit status
metadata$post$exitStatus = "complete"
# write to file
write_yaml(metadata, paste0(config$run$runDir, "/metadata/run_metadata.yml"))
