#!/bin/bash

# bookmark project folder for later reference
project_dir="$1"

# navigate to run folder
cd "$project_dir"

# prompt user to drop folders if preexisting
if [ -d "$project_dir/qc" ]; then
  echo "The existing qc folder will be deleted from the run directory ($project_dir). Proceed? (y/n)"
  read confirm
  if [[ "$confirm" == [Yy] ]]; then
    rm -r "$project_dir/qc"
    echo "Folder deleted."
  else
    echo "Process canceled."
    exit 1
  fi
fi

mkdir "$project_dir/qc"

# r01
echo "Making manifest file"
Rscript /project/source/r/00_make_manifest_file.R "$project_dir/qc" TRUE

# s01
echo "Importing Sequences as .qza file"
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path qc/manifest.csv --output-path qc/s01_16S_eDNA_Demux.qza --input-format PairedEndFastqManifestPhred33

# s02
echo "Summarizing qzv for viewing on QIIME2View"
qiime demux summarize --i-data qc/s01_16S_eDNA_Demux.qza --o-visualization qc/s02_16S_eDNA_Demux.qzv
echo "Quality control finished! Visualize quality file at https://view.qiime2.org/"