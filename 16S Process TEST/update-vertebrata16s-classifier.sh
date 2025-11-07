#!/bin/bash

# save and export config path as environmental variable
env_config_path=$config.yml
export env_config_path

echo "updating taxonomic classifier"
bash source/shell/fetch-ncbi-seq-classifier-vertebrata16s.sh

echo "writing classifier metadata"
Rscript source/r/aa_build_classifier_metadata.R "$env_config_path"