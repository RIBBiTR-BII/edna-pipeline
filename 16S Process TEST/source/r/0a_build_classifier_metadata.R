library(yaml)

# inherit env_config_path
args = commandArgs(trailingOnly = TRUE)
gene = args[1]
src = args[2]

# initiate classifier metadata
metadata = list()
metadata$classifier$name = paste0("NCBI_Vertebrata_", gene)
metadata$classifier$buildScript = src
metadata$classifier$buildEndTimestamp = format(Sys.time(), '%Y-%m-%d %H:%M:%S %z')
# write to file
write_yaml(metadata, paste0(config$taxonomy$referenceDir, "/classifier_metadata.yml"))