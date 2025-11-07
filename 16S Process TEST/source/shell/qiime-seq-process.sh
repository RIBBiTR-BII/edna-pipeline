#!/bin/bash

echo "Reading config file"
# read in config.yml (from $env_config_path) for configuration parameters
c_run_dir=$(yq e '.run.runDir' "$env_config_path")
c_taxa_dir=$(yq e '.taxonomy.classifierDir' "$env_config_path")
c_qiime_env=$(yq e '.run.qiimeEnv' "$env_config_path")
c_error_rate=$(yq e '.trimPaired.errorRate' "$env_config_path")
c_trunc_q=$(yq e '.denoisePaired.truncQ' "$env_config_path")
c_trunc_len=$(yq e '.denoisePaired.truncLen' "$env_config_path")
c_maxaccepts=$(yq e '.vsearchGlobal.maxaccepts' "$env_config_path")
c_perc_identity=$(yq e '.vsearchGlobal.percIdentity' "$env_config_path")
c_query_cov=$(yq e '.vsearchGlobal.queryCov' "$env_config_path")

# bookmark project folder for later reference
project_dir=$(pwd)

echo "Activating qiime2 environment"
# load qiime2 environment 
source activate $c_qiime_env
# visit https://docs.qiime2.org/2024.5/ to install QIIME2 & learn more about it.

# navigate to run folder
cd $c_run_dir

echo "Importing Sequences as .qza file"
# Import Sequences into qiime2 as .qza file. 
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path metadata/manifest.csv --output-path analysis/s01_16S_eDNA_Demux.qza  --input-format PairedEndFastqManifestPhred33

echo "Summarizing qzv for viewing on QIIME2View"
# Summarize into a qzv file for viewing on QIIME2View
qiime demux summarize --i-data analysis/s01_16S_eDNA_Demux.qza --o-visualization analysis/s02_16S_eDNA_Demux.qzv 

echo "Trimming primers"
# Trim primers from paired-end demux 
qiime cutadapt trim-paired --i-demultiplexed-sequences analysis/s01_16S_eDNA_Demux.qza --p-front-f ACGAGAAGACCCYRTGGARCTT --p-front-r ATCCAACATCGAGGTCGTAA --p-error-rate $c_error_rate --p-match-adapter-wildcards True --p-match-read-wildcards True --p-discard-untrimmed True --p-times 10 --o-trimmed-sequences analysis/s03_lead_primertrimmed_16S_eDNA_Demux_1.qza

# Trim primers from paired-end demux 
qiime cutadapt trim-paired --i-demultiplexed-sequences analysis/s03_lead_primertrimmed_16S_eDNA_Demux_1.qza --p-adapter-f TTACGACCTCGATGTTGGAT --p-adapter-r AAGYTCCAYRGGGTCTTCTCGT --p-error-rate 0.3 --p-match-adapter-wildcards True  --p-match-read-wildcards True  --p-discard-untrimmed False  --p-times 10  --o-trimmed-sequences analysis/s04_primertrimmed_16S_eDNA_Demux_2.qza

echo "Summarizing qzv for viewing of trimmed primers on QIIME2View"
# Summarize into a qzv file for viewing on QIIME2View
qiime demux summarize  --i-data analysis/s04_primertrimmed_16S_eDNA_Demux_2.qza  --o-visualization analysis/s05_primertrimmed_16S_eDNA_Demux.qzv 

echo "Denoising sequences with dada2"
# Denoise paired sequences with dada2
qiime dada2 denoise-paired --i-demultiplexed-seqs analysis/s04_primertrimmed_16S_eDNA_Demux_2.qza --p-trunc-q $c_trunc_q --p-trunc-len-f $c_trunc_len --p-trunc-len-r $c_trunc_len --p-n-threads 0 --output-dir analysis/s05_denoised_16S_eDNA

qiime tools extract --input-path analysis/s05_denoised_16S_eDNA/table.qza --output-path analysis/s05_denoised_16S_eDNA/sample_table

biom convert -i analysis/s05_denoised_16S_eDNA/sample_table/*/data/feature-table.biom -o analysis/s05_denoised_16S_eDNA/sample_table/feature-table.tsv --to-tsv

qiime tools extract --input-path analysis/s05_denoised_16S_eDNA/representative_sequences.qza --output-path analysis/s05_denoised_16S_eDNA/representative_sequences

# before you classify your sequences, you'll need to run the code in 'taxonomic-classification' using the rescript plug-in in qiime. 

echo "Classifying sequences"
# Classify Sequences vsearch global
qiime feature-classifier vsearch-global --i-query analysis/s05_denoised_16S_eDNA/representative_sequences.qza --i-reference-reads "$project_dir/$c_taxa_dir/Vertebrata16S_derep1_seqs_extracted.qza" --p-maxaccepts $c_maxaccepts --p-perc-identity $c_perc_identity --p-query-cov $c_query_cov --output-dir analysis/s06_classified_taxonomy_vsearch

echo "Extracting classifications"
# Extract out classifications. 
qiime tools extract  --input-path analysis/s06_classified_taxonomy_vsearch/search_results.qza --output-path analysis/s06_classified_taxonomy_vsearch/search_results

# return to project folder
echo "All finished! Proceeding to post sequencing analysis"