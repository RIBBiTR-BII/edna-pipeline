#!/bin/bash

# Activate the Qiime environment
# source activate qiime2-amplicon-2025.4

c_taxa_dir="amphibia16s_classifier"

# make sure to have rescript installed within this environment
# see details here: https://library.qiime2.org/plugins/rescript/27/

if [ -d "$c_taxa_dir" ]; then
  echo "The existing "Amphibia16s_classifier" folder will be deleted. Proceed? (y/n)"
  read confirm
  if [[ "$confirm" == [Yy] ]]; then
    rm -r "$c_taxa_dir"
    echo "Folder deleted."
  else
    echo "Process canceled."
    exit 1
  fi
fi

# make a new dirctory for the NCBI sequences
mkdir "$c_taxa_dir"
cd "$c_taxa_dir"

# download all vertebrate (txid7742) 16S sequences from NCBI
echo "Downloading Amphibian 16S sequences from NCBI"

qiime rescript get-ncbi-data \
--p-query "((txid8292[ORGN] AND (mitochondria[TITLE] OR mitochondrion[TITLE] OR mitochondrial[TITLE])) OR (txid8292[ORGN] AND (large subunit ribosomal RNA[TITLE] OR 16S rRNA[TITLE] OR 16S ribosomal RNA[TITLE] OR 16S[TITLE] OR 16S r RNA[TITLE] OR MT-RNR2[TITLE] OR MTRNR2[TITLE] OR RNR2[TITLE] OR rrnL[TITLE] OR rrn16[TITLE] OR l-rRNA[TITLE]))) AND (ddbj_embl_genbank[filter] AND mitochondrion[filter])" \
--output-dir temp

# remove sequences with 5 or more degen sequences and homopolymers longer than 12
echo "Filtering sequences"
qiime rescript cull-seqs \
    --i-sequences ./temp/sequences.qza \
    --p-num-degenerates 5 \
    --p-homopolymer-length 12 \
    --o-clean-sequences ./temp/sequences_culled.qza

# remove any replicates
echo "Removing replicates"
qiime rescript dereplicate --verbose \
  --i-sequences ./temp/sequences_culled.qza \
  --i-taxa ./temp/taxonomy.qza \
  --p-mode 'super' \
  --p-derep-prefix \
  --p-rank-handles kingdom phylum class order family genus species \
  --o-dereplicated-sequences ./temp/Amphibia16S_derep1_seqs.qza \
  --o-dereplicated-taxa ./temp/Amphibia16S_derep1_taxa.qza

# extract sequences with the primer in its entirety and remove any sequences larger or smaller than expected.
echo "Extracting sequences"
qiime feature-classifier extract-reads \
--i-sequences ./temp/Amphibia16S_derep1_seqs.qza \
--p-f-primer ACGAGAAGACCCYRTGGARCTT \
--p-r-primer ATCCAACATCGAGGTCGTAA \
--p-min-length 150 \
--p-max-length 450 \
--o-reads Amphibia16S_derep1_seqs_extracted.qza

# filter taxa to filtered sequences
qiime rescript filter-taxa \
    --i-taxonomy ./temp/Amphibia16S_derep1_taxa.qza \
    --m-ids-to-keep-file Amphibia16S_derep1_seqs_extracted.qza \
    --o-filtered-taxonomy Amphibia16S_derep1_taxa_extracted.qza

# extract taxa
qiime tools extract \
    --input-path Amphibia16S_derep1_taxa_extracted.qza \
    --output-path Amphibia16S_derep1_taxa_extracted

# extract sequences
qiime tools extract \
    --input-path Amphibia16S_derep1_seqs_extracted.qza \
    --output-path Amphibia16S_derep1_seqs_extracted

# drop temp folder
echo "deleting temporary files"
rm -r "temp"

echo "Done creating classifier"
