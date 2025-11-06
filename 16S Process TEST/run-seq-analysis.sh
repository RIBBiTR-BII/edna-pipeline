#!/bin/bash

# read in config.yml for configuration parameters
c_run_dir=$(yq e '.run.runDir' config.yml)
c_update_taxa=$(yq e '.taxonomy.updateBool' config.yml)

# prompt user to drop folders if preexisting
if [ -d $c_run_dir/metadata ]; then
  echo "The existing metadata, analysis and output folders will be deleted from the run directory ($c_run_dir). Proceed? (y/n)"
  read confirm
  if [[ "$confirm" == [Yy] ]]; then
    rm -r "$c_run_dir/metadata" "$c_run_dir/analysis" "$c_run_dir/output"
    echo "Folders deleted."
  else
    echo "Process canceled."
    exit 1
  fi
fi

echo "Creating run directories"
mkdir "$c_run_dir/metadata" "$c_run_dir/analysis" "$c_run_dir/output"

# update reference taxonomy
if [ "$c_update_taxa" = "true" ]; then
  echo "Downloading and curating reference taxa list"
  bash update-vertebrata16s-classifier.sh 
fi

echo "Making manifest file"
Rscript source/r/00_make_manifest_file.R

echo "Initiating run metadata"
Rscript source/r/01_init_run_metadata.R

pwd

echo "performing sequence processing"
bash source/shell/qiime-seq-process.sh

pwd

echo "Finalizing run metadata"
Rscript source/r/02_final_run_metadata.R

# echo "performing sequence analysis"
# Rscript source/r/03_consolidate_outputs.R

echo "Sequence analysis complete! See output folder for relevant outputs"