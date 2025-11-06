library(yaml)
library(lubridate)

# read in config file, use as basis for metadata
metadata = read_yaml("config.yml")
# scrape classifier metadata
classifier_md = read_yaml(paste0(metadata$taxonomy$referenceDir, "/classifier_metadata.yml"))
metadata$classifier = classifier_md$classifier
# scrape manifest file
manifest = read.csv(paste0(metadata$run$runDir, "/metadata/manifest.csv"))
metadata$sequences$sampleCount = length(unique(manifest$sample.id))
metadata$sequences$sampleIDFirst = head(manifest$sample.id, 1)
metadata$sequences$sampleIDLast = tail(manifest$sample.id, 1)
# add start and end time
metadata$runtime$startTimestamp = format(now(), '%Y-%m-%d %H:%M:%S')
metadata$runtime$endTimestamp = NA
#write to file
write_yaml(metadata, paste0(config$run$runDir, "/metadata/run_metadata.yml"))
