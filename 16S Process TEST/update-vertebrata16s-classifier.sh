#!/bin/bash

echo "updating taxonomic classifier"
. source/shell/fetch-ncbi-seq-classifier-vertebrata16s.sh

echo "writing classifier metadata"
Rscript source/r/aa_buld_classifier_metadata.R