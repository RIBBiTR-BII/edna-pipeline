#!/bin/bash

echo "updating taxonomic classifier"
bash source/shell/fetch-ncbi-seq-classifier-vertebrata16s.sh

echo "writing classifier metadata"
Rscript source/r/aa_build_classifier_metadata.R