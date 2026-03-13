#!/bin/bash

echo "Reading config file"
# read in config.yml (from $env_config_path) for configuration parameters
c_run_dir=$(yq e '.run.runDir' "$env_config_path")
c_taxa_dir=$(yq e '.taxonomy.classifierDir' "$env_config_path")
c_qiime_env=$(yq e '.run.qiimeEnv' "$env_config_path")
c_error_rate=$(yq e '.trimPaired.errorRate' "$env_config_path")
c_trunc_q=$(yq e '.denoisePaired.truncQ' "$env_config_path")
c_trunc_len=$(yq e '.denoisePaired.truncLen' "$env_config_path")
c_b_maxaccepts=$(yq e '.consensusBlast.maxaccepts' "$env_config_path")
c_b_perc_identity=$(yq e '.consensusBlast.percIdentity' "$env_config_path")
c_b_query_cov=$(yq e '.consensusBlast.queryCov' "$env_config_path")
c_v_maxaccepts=$(yq e '.vsearchGlobal.maxaccepts' "$env_config_path")
c_v_perc_identity=$(yq e '.vsearchGlobal.percIdentity' "$env_config_path")
c_v_query_cov=$(yq e '.vsearchGlobal.queryCov' "$env_config_path")

# bookmark project folder for later reference
project_dir=$(pwd)\

# load qiime2 environment 
echo "Activating qiime2 environment"
source activate $c_qiime_env
# visit https://docs.qiime2.org/2024.5/ to install QIIME2 & learn more about it.

# navigate to run folder
cd $c_run_dir

# s01
# Import Sequences into qiime2 as .qza file. 
echo "Importing Sequences as .qza file"
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path output/metadata/manifest.csv --output-path analysis/s01_12S_eDNA_Demux.qza  --input-format PairedEndFastqManifestPhred33

# s02
# Summarize into a qzv file for viewing on QIIME2View
echo "Summarizing qzv for viewing on QIIME2View"
qiime demux summarize --i-data analysis/s01_12S_eDNA_Demux.qza --o-visualization analysis/s02_12S_eDNA_Demux.qzv 

# s03
# Trim primers from paired-end demux (lead)
echo "Trimming primers"
qiime cutadapt trim-paired --i-demultiplexed-sequences analysis/s01_12S_eDNA_Demux.qza --p-front-f ACACCGCCCGTCACCCT --p-front-r GTAYACTTACCATGTTACGACT --p-error-rate $c_error_rate --p-match-adapter-wildcards True --p-match-read-wildcards True --p-discard-untrimmed True --p-times 10 --o-trimmed-sequences analysis/s03_lead_primertrimmed_12S_eDNA_Demux.qza

# s04
# Trim primers from paired-end demux (lag)
qiime cutadapt trim-paired --i-demultiplexed-sequences analysis/s03_lead_primertrimmed_12S_eDNA_Demux.qza --p-adapter-f AAGTCGTAACATGGTAAGTRTA --p-adapter-r AGGGTGACGGGCGGTG --p-error-rate 0.3 --p-match-adapter-wildcards True  --p-match-read-wildcards True  --p-discard-untrimmed False  --p-times 10  --o-trimmed-sequences analysis/s04_leadlag_primertrimmed_12S_eDNA_Demux.qza

# s05
# Summarize into a qzv file for viewing on QIIME2View
echo "Summarizing qzv for viewing of trimmed primers on QIIME2View"
qiime demux summarize  --i-data analysis/s04_leadlag_primertrimmed_12S_eDNA_Demux.qza  --o-visualization analysis/s05_primertrimmed_12S_eDNA_Demux.qzv 

# s06
# Denoise paired sequences with dada2
echo "Denoising sequences with dada2"
qiime dada2 denoise-paired --i-demultiplexed-seqs analysis/s04_leadlag_primertrimmed_12S_eDNA_Demux.qza --p-trunc-q $c_trunc_q --p-trunc-len-f $c_trunc_len --p-trunc-len-r $c_trunc_len --p-n-threads 0 --output-dir analysis/s06_denoised_12S_eDNA

qiime tools extract --input-path analysis/s06_denoised_12S_eDNA/table.qza --output-path analysis/s06_denoised_12S_eDNA/sample_table

biom convert -i analysis/s06_denoised_12S_eDNA/sample_table/*/data/feature-table.biom -o analysis/s06_denoised_12S_eDNA/sample_table/feature-table.tsv --to-tsv

qiime tools extract --input-path analysis/s06_denoised_12S_eDNA/representative_sequences.qza --output-path analysis/s06_denoised_12S_eDNA/representative_sequences


echo "Classifying sequences"

# s07-a
# Classify Sequences consensus blast
qiime feature-classifier classify-consensus-blast --i-query analysis/s06_denoised_12S_eDNA/representative_sequences.qza --i-reference-reads "$project_dir/$c_taxa_dir/Vertebrata12S_derep1_seqs_extracted.qza" --i-reference-taxonomy "$project_dir/$c_taxa_dir/Vertebrata12S_derep1_taxa_extracted.qza" --p-maxaccepts $c_b_maxaccepts --p-perc-identity $c_b_perc_identity --p-query-cov $c_b_query_cov --output-dir analysis/s07_classified_taxonomy_blast
# Extract out classifications.
qiime tools extract  --input-path analysis/s07_classified_taxonomy_blast/classification.qza --output-path analysis/s07_classified_taxonomy_blast/classification
qiime tools extract  --input-path analysis/s07_classified_taxonomy_blast/search_results.qza --output-path analysis/s07_classified_taxonomy_blast/search_results

# s07-b
# Classify Sequences vsearch global
qiime feature-classifier vsearch-global --i-query analysis/s06_denoised_12S_eDNA/representative_sequences.qza --i-reference-reads "$project_dir/$c_taxa_dir/Vertebrata12S_derep1_seqs_extracted.qza" --p-maxaccepts $c_v_maxaccepts --p-perc-identity $c_v_perc_identity --p-query-cov $c_v_query_cov --output-dir analysis/s07_classified_taxonomy_vsearch
# Extract out classifications. 
qiime tools extract  --input-path analysis/s07_classified_taxonomy_vsearch/search_results.qza --output-path analysis/s07_classified_taxonomy_vsearch/search_results

# return to project folder
echo "All finished! Proceeding to post sequencing analysis"