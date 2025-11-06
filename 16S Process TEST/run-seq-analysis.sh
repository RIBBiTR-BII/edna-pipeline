#!/bin/bash

# read in config.yml for configuration parameters
c_run_dir=$(yq e '.run.runDir' config.yml)
c_update_taxa=$(yq e '.taxonomy.updateBool' config.yml)

# prompt user to drop folders if preexisting
if [ -d $c_run_dir/metadata ]; then
  echo "The existing "metadata", "analysis" and "output" folders will be deleted. Proceed? (y/n)"
  read confirm
  if [[ "$confirm" == [Yy] ]]; then
    rm -r $c_run_dir/metadata $c_run_dir/analysis $c_run_dir/output
    echo "Folders deleted."
  else
    echo "Process canceled."
    exit 1
  fi
fi

echo "creating necessary directories"
mkdir $c_run_dir/metadata $c_run_dir/analysis $c_run_dir/output

# update reference taxonomy
if [ $c_update_taxa ]; then
  echo "downloading and curating reference taxa list"
  . source/shell/DownloadSeqs_And_TrainClassifier_Vertebrata16S.sh 
fi

echo "making manifest file"
Rscript source/r/00_make_manifest_file.R

echo "initiating run metadata"
Rscript source/r/01_init_run_metadata.R

echo "performing sequence processing"
. source/shell/qiime-seq-process.sh

echo "finalizing run metadata"
Rscript source/r/02_final_run_metadata.R

echo "performing sequence analysis"
Rscript source/r/02_consolidate_outputs.R

echo "sequence analysis complete.  please see output folder for relevant output"