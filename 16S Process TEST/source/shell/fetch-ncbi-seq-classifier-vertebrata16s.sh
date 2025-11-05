#!/bin/bash

# Activate the Qiime environment
source activate qiime2-amplicon-2025.4

# make sure to have rescript installed within this environment
# see details here: https://library.qiime2.org/plugins/rescript/27/

if [ -d "vertebrata16s-classifier" ]; then
  echo "The existing "vertebrata16s-classifier" folder will be deleted. Proceed? (y/n)"
  read confirm
  if [[ "$confirm" == [Yy] ]]; then
    rm -r "vertebrata16s-classifier"
    echo "Folder deleted."
  else
    echo "Process canceled."
    exit 1
  fi
fi

# make a new dirctory for the NCBI sequences
mkdir vertebrata16s-classifier
cd vertebrata16s-classifier

# download all vertebrate (txid7742) 16S sequences from NCBI
echo "downloading vertibrate 16S sequences from NCBI"
qiime rescript get-ncbi-data \
--p-query "(txid7742[ORGN] AND (mitochondria[TITLE] OR mitochondrion[TITLE] OR mitochondrial[TITLE])) OR (txid7742[ORGN] AND (large subunit ribosomal RNA[TITLE] OR 16S rRNA[TITLE] OR 16S ribosomal RNA[TITLE] OR 16S[TITLE] OR 16S r RNA[TITLE] OR MT-RNR2[TITLE] OR MTRNR2[TITLE] OR RNR2[TITLE])) AND (biomol_genomic[PROP] AND ddbj_embl_genbank[filter] AND mitochondrion[filter])" \
--output-dir temp

# remove sequences with 5 or more degen sequences and homopolymers longer than 12
echo "filtering sequences"
qiime rescript cull-seqs \
    --i-sequences ./temp/sequences.qza \
    --p-num-degenerates 5 \
    --p-homopolymer-length 12 \
    --o-clean-sequences ./temp/Vertebrata16S_ambi_hpoly_filtd_seqs.qza

# remove any replicates
echo "removing replicates"
qiime rescript dereplicate --verbose \
  --i-sequences ./temp/sequences.qza \
  --i-taxa ./temp/taxonomy.qza \
  --p-mode 'super' \
  --p-derep-prefix \
  --p-rank-handles ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'] \
  --o-dereplicated-sequences ./temp/Vertebrata16S_derep1_seqs.qza \
  --o-dereplicated-taxa ./temp/Vertebrata16S_derep1_taxa.qza

# extract sequences with the primer in its entirety and remove any sequences larger or smaller than expected.
echo "extracting sequences"
qiime feature-classifier extract-reads \
--i-sequences ./temp/Vertebrata16S_derep1_seqs.qza \
--p-f-primer ACGAGAAGACCCYRTGGARCTT \
--p-r-primer ATCCAACATCGAGGTCGTAA \
--p-min-length 150 \
--p-max-length 450 \
--o-reads Vertebrata16S_derep1_seqs_extracted.qza

qiime rescript filter-taxa \
    --i-taxonomy ./temp/Vertebrata16S_derep1_taxa.qza \
    --m-ids-to-keep-file Vertebrata16S_derep1_seqs_extracted.qza \
    --o-filtered-taxonomy Vertebrata16S_derep1_taxa_extracted.qza

# # drop temp folder
# echo "deleting temporary files"
# rm -r "temp"

echo "done creating classifier"
