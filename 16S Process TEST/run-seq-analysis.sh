#!/bin/bash

# read in config.yml for configuration parameters
c_run_dir=$(yq e '.run-directory' config.yml)

echo "clearing previous folders that will be created by this pipeline"
rm -r $c_run_dir/output $c_run_dir/metadata $c_run_dir/analysis


echo "creating necessary directories"
mkdir $c_run_dir/output $c_run_dir/metadata $c_run_dir/analysis

#echo "downloading and curating reference taxa list"
#. source/shell/DownloadSeqs_And_TrainClassifier_Vertebrata16S.sh 

echo "making manifest file"
Rscript source/r/00_make_manifest_file.R

echo "performing sequence processing"
. source/shell/qiime-seq-process.sh

echo "performing sequence analysis"
Rscript source/r/01_consolidate_outputs.R

echo "sequence analysis complete.  please see output folder for relevant output"