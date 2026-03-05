#!/bin/bash

gene="12S"
src="./source/shell/fetch-ncbi-seq-classifier-vertebrata12s.sh"

echo "updating taxonomic classifier"
bash source/shell/fetch-ncbi-seq-classifier-vertebrata12s.sh

echo "writing classifier metadata"
Rscript source/r/0a_build_classifier_metadata.R "$gene" "$src"