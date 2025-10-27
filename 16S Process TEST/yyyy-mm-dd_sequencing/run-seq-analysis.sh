#!/bin/bash

echo "clearing previous folders that will be created by this pipeline"
rm -r output metadata analysis


echo "creating necessary directories"
mkdir output metadata analysis

#echo "downloading and curating reference taxa list"
#sh source/shell/DownloadSeqs_And_TrainClassifier_Vertebrata16S.sh 

echo "making manifest file"
Rscript source/r/00_make_manifest_file.R

echo "performing sequence processing"
. source/shell/qiime-seq-process.sh

echo "performing sequence analysis"
Rscript source/r/01_making_relevant_files.R

echo "sequence analysis complete.  please see output folder for relevant output"