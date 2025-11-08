#!/bin/bash

echo "activating qiime2 environment"
# load qiime2 environment 
source activate qiime2-amplicon-2024.5

#go to https://docs.qiime2.org/2024.5/ to install QIIME2 & learn more about it. 

#use ls navigate to Example_eDNA_Pipeline and cd to set it as working directory

# make sure you have folders named 'metadata', 'output', 'sequences', and 'taxonomic-classification' in this folder 

echo "Importing Sequences as .qza file"
# Import Sequences into qiime2 as .qza file. 
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path metadata/manifest.csv --output-path analysis/16S_eDNA_Demux.qza  --input-format PairedEndFastqManifestPhred33

echo "Summarizing qzv for viewing on QIIME2View"
# Summarize into a qzv file for viewing on QIIME2View
qiime demux summarize --i-data analysis/16S_eDNA_Demux.qza --o-visualization analysis/16S_eDNA_Demux.qzv 

echo "Trimming primers"
# Trim primers from paired-end demux 
qiime cutadapt trim-paired --i-demultiplexed-sequences analysis/16S_eDNA_Demux.qza --p-front-f ACGAGAAGACCCYRTGGARCTT --p-front-r ATCCAACATCGAGGTCGTAA --p-error-rate 0.3 --p-match-adapter-wildcards True --p-match-read-wildcards True --p-discard-untrimmed True --p-times 10 --o-trimmed-sequences analysis/lead_primertrimmed_16S_eDNA_Demux.qza

# Trim primers from paired-end demux 
qiime cutadapt trim-paired --i-demultiplexed-sequences analysis/lead_primertrimmed_16S_eDNA_Demux.qza --p-adapter-f TTACGACCTCGATGTTGGAT --p-adapter-r AAGYTCCAYRGGGTCTTCTCGT --p-error-rate 0.3 --p-match-adapter-wildcards True  --p-match-read-wildcards True  --p-discard-untrimmed False  --p-times 10  --o-trimmed-sequences analysis/primertrimmed_16S_eDNA_Demux.qza

echo "Summarizing qzv for viewing of trimmed primers on QIIME2View"
# Summarize into a qzv file for viewing on QIIME2View
qiime demux summarize  --i-data analysis/primertrimmed_16S_eDNA_Demux.qza  --o-visualization analysis/primertrimmed_16S_eDNA_Demux.qzv 

echo "Denoising sequences with dada2"
# Denoise paired sequences with dada2
qiime dada2 denoise-paired --i-demultiplexed-seqs analysis/primertrimmed_16S_eDNA_Demux.qza --p-trunc-q 2 --p-trunc-len-f 0 --p-trunc-len-r 0 --p-n-threads 0 --output-dir analysis/denoised_16S_eDNA

qiime tools extract --input-path analysis/denoised_16S_eDNA/table.qza --output-path analysis/denoised_16S_eDNA/sample_table

biom convert -i analysis/denoised_16S_eDNA/sample_table/*/data/feature-table.biom -o analysis/denoised_16S_eDNA/sample_table/feature-table.tsv --to-tsv

qiime tools extract --input-path analysis/denoised_16S_eDNA/representative_sequences.qza --output-path analysis/denoised_16S_eDNA/representative_sequences

# before you classify your sequences, you'll need to run the code in 'taxonomic-classification' using the rescript plug-in in qiime. 

echo "Classifying sequences"
# Classify Sequences
qiime feature-classifier classify-consensus-blast --i-query analysis/denoised_16S_eDNA/representative_sequences.qza --i-reference-reads taxonomic-classification/Vertebrata16S_derep1_seqs.qza --i-reference-taxonomy taxonomic-classification/Vertebrata16S_derep1_taxa.qza --p-maxaccepts 10 --p-perc-identity 0.98 --p-query-cov 0.95 --output-dir analysis/denoised_16S_eDNA/classified_taxonomy

echo "Extracting classifications"
# Extract out classifications. 
qiime tools extract  --input-path analysis/denoised_16S_eDNA/classified_taxonomy/classification.qza --output-path analysis/denoised_16S_eDNA/classified_taxonomy/classified_taxonomy

echo "All finished! Proceeding to post sequencing analysis"