#!/bin/bash

echo "reading config file"
# read in config.yml for configuration parameters
c_run_dir=$(yq e '.run-directory' config.yml)
c_error_rate=$(yq e '.trim-paired.p-error-rate' config.yml)
c_trunc_q=$(yq e '.denoise-paired.p-trunc-q' config.yml)
c_trunc_len=$(yq e '.denoise-paired.p-trunc-len' config.yml)
c_maxaccepts=$(yq e '.clasify-consensus-blast.p-maxaccepts' config.yml)
c_perc_identity=$(yq e '.clasify-consensus-blast.p-perc-identity' config.yml)
c_query_cov=$(yq e '.clasify-consensus-blast.p-query-cov' config.yml)

echo "activating qiime2 environment"
# load qiime2 environment 
source activate qiime2-amplicon-2025.4

#go to https://docs.qiime2.org/2024.5/ to install QIIME2 & learn more about it. 

# use ls navigate to Example_eDNA_Pipeline and cd to set it as working directory
project_dir=$(pwd)
cd $c_run_dir

# make sure you have folders named 'metadata', 'output', 'sequences', and 'taxonomic-classification' in this folder 

echo "Importing Sequences as .qza file"
# Import Sequences into qiime2 as .qza file. 
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path metadata/manifest.csv --output-path analysis/16S_eDNA_Demux.qza  --input-format PairedEndFastqManifestPhred33

echo "Summarizing qzv for viewing on QIIME2View"
# Summarize into a qzv file for viewing on QIIME2View
qiime demux summarize --i-data analysis/16S_eDNA_Demux.qza --o-visualization analysis/16S_eDNA_Demux.qzv 

echo "Trimming primers"
# Trim primers from paired-end demux 
qiime cutadapt trim-paired --i-demultiplexed-sequences analysis/16S_eDNA_Demux.qza --p-front-f ACGAGAAGACCCYRTGGARCTT --p-front-r ATCCAACATCGAGGTCGTAA --p-error-rate $c_error_rate --p-match-adapter-wildcards True --p-match-read-wildcards True --p-discard-untrimmed True --p-times 10 --o-trimmed-sequences analysis/lead_primertrimmed_16S_eDNA_Demux.qza

# Trim primers from paired-end demux 
qiime cutadapt trim-paired --i-demultiplexed-sequences analysis/lead_primertrimmed_16S_eDNA_Demux.qza --p-adapter-f TTACGACCTCGATGTTGGAT --p-adapter-r AAGYTCCAYRGGGTCTTCTCGT --p-error-rate 0.3 --p-match-adapter-wildcards True  --p-match-read-wildcards True  --p-discard-untrimmed False  --p-times 10  --o-trimmed-sequences analysis/primertrimmed_16S_eDNA_Demux.qza

echo "Summarizing qzv for viewing of trimmed primers on QIIME2View"
# Summarize into a qzv file for viewing on QIIME2View
qiime demux summarize  --i-data analysis/primertrimmed_16S_eDNA_Demux.qza  --o-visualization analysis/primertrimmed_16S_eDNA_Demux.qzv 

echo "Denoising sequences with dada2"
# Denoise paired sequences with dada2
qiime dada2 denoise-paired --i-demultiplexed-seqs analysis/primertrimmed_16S_eDNA_Demux.qza --p-trunc-q $c_trunc_q --p-trunc-len-f $c_trunc_len --p-trunc-len-r $c_trunc_len --p-n-threads 0 --output-dir analysis/denoised_16S_eDNA

qiime tools extract --input-path analysis/denoised_16S_eDNA/table.qza --output-path analysis/denoised_16S_eDNA/sample_table

biom convert -i analysis/denoised_16S_eDNA/sample_table/*/data/feature-table.biom -o analysis/denoised_16S_eDNA/sample_table/feature-table.tsv --to-tsv

qiime tools extract --input-path analysis/denoised_16S_eDNA/representative_sequences.qza --output-path analysis/denoised_16S_eDNA/representative_sequences

# before you classify your sequences, you'll need to run the code in 'taxonomic-classification' using the rescript plug-in in qiime. 

echo "Classifying sequences"
# Classify Sequences
# consensus BLAST
qiime feature-classifier classify-consensus-blast --i-query analysis/denoised_16S_eDNA/representative_sequences.qza --i-reference-reads $project_dir/taxonomic-classification/Vertebrata16S_derep1_seqs.qza --i-reference-taxonomy $project_dir/taxonomic-classification/Vertebrata16S_derep1_taxa.qza --p-maxaccepts $c_maxaccepts --p-perc-identity $c_perc_identity --p-query-cov $c_query_cov --output-dir analysis/denoised_16S_eDNA/classified_taxonomy

# vsearch global
# qiime feature-classifier vsearch-global --i-query analysis/denoised_16S_eDNA/representative_sequences.qza --i-reference-reads taxonomic-classification/Vertebrata16S_derep1_seqs.qza --p-maxaccepts 20 --p-perc-identity 0.90 --p-query-cov 0.70 --output-dir analysis/denoised_16S_eDNA/classified_taxonomy_vsearch

echo "Extracting classifications"
# Extract out classifications. 
qiime tools extract  --input-path $c_run_dir/analysis/denoised_16S_eDNA/classified_taxonomy/classification.qza --output-path $c_run_dir/analysis/denoised_16S_eDNA/classified_taxonomy/classified_taxonomy

echo "All finished! Proceeding to post sequencing analysis"