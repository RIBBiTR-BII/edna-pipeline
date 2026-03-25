#!/bin/bash

gene="16S"
src="./source/shell/fetch-ncbi-seq-classifier-amphibia16s.sh"

echo "updating taxonomic classifier"
bash source/shell/fetch-ncbi-seq-classifier-amphibia16s.sh

echo "writing classifier metadata"
Rscript source/r/0a_build_classifier_metadata.R "$gene" "$src"