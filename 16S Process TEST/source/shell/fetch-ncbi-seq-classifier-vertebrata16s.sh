#!/bin/bash

# read in config.yml (from $env_config_path) for configuration parameters
c_taxa_dir=$(yq e '.taxonomy.classifierDir' "$env_config_path")
c_qiime_env=$(yq e '.run.qiimeEnv' "$env_config_path")

# Activate the Qiime environment
source activate $c_qiime_env

# make sure to have rescript installed within this environment
# see details here: https://library.qiime2.org/plugins/rescript/27/

if [ -d $c_taxa_dir ]; then
  echo "The existing "vertebrata16s-classifier" folder will be deleted. Proceed? (y/n)"
  read confirm
  if [[ "$confirm" == [Yy] ]]; then
    rm -r $c_taxa_dir
    echo "Folder deleted."
  else
    echo "Process canceled."
    exit 1
  fi
fi

# make a new dirctory for the NCBI sequences
mkdir $c_taxa_dir
cd $c_taxa_dir

# download all vertebrate (txid7742) 16S sequences from NCBI
echo "Downloading vertibrate 16S sequences from NCBI"
qiime rescript get-ncbi-data \
--p-query "(txid7742[ORGN] AND (mitochondria[TITLE] OR mitochondrion[TITLE] OR mitochondrial[TITLE])) OR (txid7742[ORGN] AND (large subunit ribosomal RNA[TITLE] OR 16S rRNA[TITLE] OR 16S ribosomal RNA[TITLE] OR 16S[TITLE] OR 16S r RNA[TITLE] OR MT-RNR2[TITLE] OR MTRNR2[TITLE] OR RNR2[TITLE])) AND (biomol_genomic[PROP] AND ddbj_embl_genbank[filter] AND mitochondrion[filter])" \
--output-dir temp

# remove sequences with 5 or more degen sequences and homopolymers longer than 12
echo "Filtering sequences"
qiime rescript cull-seqs \
    --i-sequences ./temp/sequences.qza \
    --p-num-degenerates 5 \
    --p-homopolymer-length 12 \
    --o-clean-sequences ./temp/Vertebrata16S_ambi_hpoly_filtd_seqs.qza

# remove any replicates
echo "Removing replicates"
qiime rescript dereplicate --verbose \
  --i-sequences ./temp/sequences.qza \
  --i-taxa ./temp/taxonomy.qza \
  --p-mode 'super' \
  --p-derep-prefix \
  --p-rank-handles kingdom phylum class order family genus species \
  --o-dereplicated-sequences ./temp/Vertebrata16S_derep1_seqs.qza \
  --o-dereplicated-taxa ./temp/Vertebrata16S_derep1_taxa.qza

# extract sequences with the primer in its entirety and remove any sequences larger or smaller than expected.
echo "Extracting sequences"
qiime feature-classifier extract-reads \
--i-sequences ./temp/Vertebrata16S_derep1_seqs.qza \
--p-f-primer ACGAGAAGACCCYRTGGARCTT \
--p-r-primer ATCCAACATCGAGGTCGTAA \
--p-min-length 150 \
--p-max-length 450 \
--o-reads Vertebrata16S_derep1_seqs_extracted.qza

# filter taxa to filtered sequences
qiime rescript filter-taxa \
    --i-taxonomy ./temp/Vertebrata16S_derep1_taxa.qza \
    --m-ids-to-keep-file Vertebrata16S_derep1_seqs_extracted.qza \
    --o-filtered-taxonomy Vertebrata16S_derep1_taxa_extracted.qza

# extract taxa
qiime tools extract \
    --input-path Vertebrata16S_derep1_taxa_extracted.qza \
    --output-path Vertebrata16S_derep1_taxa_extracted

# # drop temp folder
# echo "deleting temporary files"
# rm -r "temp"

echo "Done creating classifier"
