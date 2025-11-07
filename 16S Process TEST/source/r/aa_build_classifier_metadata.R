library(yaml)

# inherit env_config_path
args = commandArgs(trailingOnly = TRUE)
env_config_path = args[1]

# read in config
config = read_yaml(env_config_path)

# initiate classifier metadata
metadata = list()
metadata$classifier$name = "NCBI_Vertebrata16S"
metadata$classifier$buildScript = "./source/shell/fetch-ncbi-seq-classifier-vertebrata16s.sh"
metadata$classifier$buildEndTimestamp = format(Sys.time(), '%Y-%m-%d %H:%M:%S %z')
# write to file
write_yaml(metadata, paste0(config$taxonomy$referenceDir, "/classifier_metadata.yml"))