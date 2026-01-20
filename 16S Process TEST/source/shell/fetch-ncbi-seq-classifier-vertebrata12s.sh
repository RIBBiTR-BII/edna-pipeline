#!/bin/bash

# Activate the Qiime environment
source activate qiime2-amplicon-2025.4

c_taxa_dir="vertebrata12s_classifier"

# make sure to have rescript installed within this environment
# see details here: https://library.qiime2.org/plugins/rescript/27/

if [ -d "$c_taxa_dir" ]; then
  echo "The existing "vertebrata12s-classifier" folder will be deleted. Proceed? (y/n)"
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


# download all vertebrate (txid7777) 12S sequences from NCBI
# echo "Downloading vertibrate 12S sequences from NCBI"
# qiime rescript get-ncbi-data \
#   --p-query "(txid7777[ORGN] AND (mitochondria[TITLE] OR mitochondrion[TITLE] OR mitochondrial[TITLE])) OR (txid7777[ORGN] AND (small subunit ribosomal RNA[TITLE] OR 12S rRNA[TITLE] OR 12S ribosomal RNA[TITLE] OR 12S[TITLE] OR 12S r RNA[TITLE] OR MT-RNR1[TITLE] OR MTRNR1[TITLE] OR RNR1[TITLE])) AND (biomol_genomic[PROP] AND ddbj_embl_genbank[filter] AND mitochondrion[filter])" \
#   --output-dir temp

# echo "Downloading vertibrate 12S sequences from NCBI"
# qiime rescript get-ncbi-data \
#   --p-query '12S[title] OR 12S[All] OR 12S ribosomal RNA[title] OR 12S[Gene] AND txid7742[ORGN]' \
#   --output-dir temp

echo "Downloading vertibrate 12S sequences from NCBI"
qiime rescript get-ncbi-data \
  --p-query " txid7742[ORGN] AND (12S OR 12S ribosomal RNA OR 12S rRNA) AND (mitochondrion[Filter] OR plastid[Filter]) NOT environmental sample[Title] NOT environmental samples[Title] NOT environmental[Title] NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title]" \
  --output-dir temp


# remove sequences with 5 or more degen sequences and homopolymers longer than 12
echo "Filtering sequences"
qiime rescript cull-seqs \
    --i-sequences ./temp/sequences.qza \
    --p-num-degenerates 5 \
    --p-homopolymer-length 12 \
    --o-clean-sequences ./temp/Vertebrata12S_ambi_hpoly_filtd_seqs.qza

# remove any replicates
echo "Removing replicates"
qiime rescript dereplicate --verbose \
  --i-sequences ./temp/sequences.qza \
  --i-taxa ./temp/taxonomy.qza \
  --p-mode 'super' \
  --p-derep-prefix \
  --p-rank-handles kingdom phylum class order family genus species \
  --o-dereplicated-sequences ./temp/Vertebrata12S_derep1_seqs.qza \
  --o-dereplicated-taxa ./temp/Vertebrata12S_derep1_taxa.qza

# extract sequences with the primer in its entirety and remove any sequences larger or smaller than expected.
echo "Extracting sequences"
qiime feature-classifier extract-reads \
--i-sequences ./temp/Vertebrata12S_derep1_seqs.qza \
--p-f-primer ACACCGCCCGTCACCCT \
--p-r-primer GTAYACTTACCATGTTACGACTT \
--p-min-length 30 \
--p-max-length 80 \
--o-reads Vertebrata12S_derep1_seqs_extracted.qza

# filter taxa to filtered sequences
qiime rescript filter-taxa \
    --i-taxonomy ./temp/Vertebrata12S_derep1_taxa.qza \
    --m-ids-to-keep-file Vertebrata12S_derep1_seqs_extracted.qza \
    --o-filtered-taxonomy Vertebrata12S_derep1_taxa_extracted.qza

# extract taxa
qiime tools extract \
    --input-path Vertebrata12S_derep1_taxa_extracted.qza \
    --output-path Vertebrata12S_derep1_taxa_extracted

# drop temp folder
echo "deleting temporary files"
rm -r "temp"

echo "Done creating classifier"
