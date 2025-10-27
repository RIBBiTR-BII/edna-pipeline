# Activate the Qiime2-2022.2 environment
source activate qiime2-2022.2

# make sure to have rescript installed within this environment
# see details here: https://library.qiime2.org/plugins/rescript/27/

cd analysis

# make a new dirctory for the NCBI sequences
mkdir 16S_classifier

# change to NCBI sequences directory
cd 16S_classifier

# Download all vertebrate (txid7742) 16S sequences from NCBI
qiime rescript get-ncbi-data \
--p-query "(txid7742[ORGN] AND (mitochondria[TITLE] OR mitochondrion[TITLE] OR mitochondrial[TITLE])) OR (txid7742[ORGN] AND (large subunit ribosomal RNA[TITLE] OR 16S rRNA[TITLE] OR 16S ribosomal RNA[TITLE] OR 16S[TITLE] OR 16S r RNA[TITLE] OR MT-RNR2[TITLE] OR MTRNR2[TITLE] OR RNR2[TITLE])) AND (biomol_genomic[PROP] AND ddbj_embl_genbank[filter] AND mitochondrion[filter])" \
--output-dir Vertebrata16S

#remove sequences with 5 or more degen sequences and homopolymers longer than 12
qiime rescript cull-seqs \
    --i-sequences ./Vertebrata16S/sequences.qza \
    --p-num-degenerates 5 \
    --p-homopolymer-length 12 \
    --o-clean-sequences Vertebrata16S_ambi_hpoly_filtd_seqs.qza

# remove any replicates
qiime rescript dereplicate --verbose \
  --i-sequences ./Vertebrata16S/sequences.qza \
  --i-taxa ./Vertebrata16S/taxonomy.qza \
  --p-mode 'super' \
  --p-derep-prefix \
  --p-rank-handles 'silva' \
  --o-dereplicated-sequences Vertebrata16S_derep1_seqs.qza \
  --o-dereplicated-taxa Vertebrata16S_derep1_taxa.qza

# extract sequences with the primer in its entirety and remove any sequences larger or smaller than expected. 
qiime feature-classifier extract-reads \
--i-sequences Vertebrata16S_derep1_seqs.qza \
--p-f-primer ACGAGAAGACCCYRTGGARCTT \
--p-r-primer ATCCAACATCGAGGTCGTAA \
--p-min-length 150 \
--p-max-length 450 \
--o-reads Vertebrata16S_derep1_seqs_extracted.qza

# extract reads to check
qiime tools extract \
--input-path Vertebrata16S_derep1_seqs_extracted.qza \
--output-path Vertebrata16S_derep1_seqs_extracted

qiime tools extract \
--input-path Vertebrata16S_derep1_taxa.qza \
--output-path Vertebrata16S_derep1_taxa 


