Title: Example_eDNA_Pipeline
Created by: Brandon Hoenig (brandonhoenig@gmail.com)


This folder will walk you through all of the steps required for analyzing eDNA in the RZ Lab. Do not use this file for bioinformatics - use Sequence_Processing_Protocol.rtf or something like it.  


NOTE: The specifics of the bioinformatics steps are *_still in development_* as options for bioninformatic processing are quite subjective.  Therefore, you may find that other bioinformatic pipelines in the data that I have passed off may not match 100% with those found here.  However, it should be noted that each of these are valid methods for processing the data and the resulting data is to be trusted.  

I do encourage the user of this file to critically observe each decision that I made to ensure that it is sound and to make updates as necessary. 

Now the guide:

First, make sure you have QIIME2 installed on your machine
  - https://docs.qiime2.org/2024.5/install/
  
You'll also want to have worked your way through the Moving Pictures Tutorial to learn more about this process. 
  - https://docs.qiime2.org/2024.5/tutorials/moving-pictures/
  
You'll also need to have R installed and have these libraries on your machine:
  qiime2R
  tidyverse
  janitor
  phylotools
  readxl
  openxlsx
  
Now to the step-by-step guide!
  
**Preparation** 
  
1: Make an all-encompassing folder. Ours is called 'Example_eDNA_Pipeline' but it can be anything you'd like really.  I like 'Month_Year_Sequencing' (e.g., August_2024_Sequencing)

2: Within this all-encompassin folder, create folders called 'metadata', 'output', 'r-code', 'r-output', 'sequences', and 'taxonomic-classification'

2A: Feel free to use the 'taxonomic-classification' folder found here as it already has the files needed. Otherwise, you'll have to run the .sh file within this folder to make your reference files for taxonomic classification. 

3: Open your terminal if on mac or your virtual box terminal in linux if you're on windows. (see QIIME2 installation guide above). 

**Sequence Processing** 

NOTE: Now that your terminal (assuming mac for the rest of this guide) is open, we will go step by step with each function in the 'Sequence_Processing_Protocol.rtf' 

NOTE: You can copy all of the text in this file into your terminal and it will run all of the steps in order.  For the first time, I'd go one by one though. 
  
1. First, you'll want to make sure you are in a bash or zsh environment. 

2. Next, navigate to the 'all-encompassing folder' by using ls (list files) and cd (change directory.). **YOUR WORKING DIRECTORY MUST BE THE ALL-ENCOMPASSING FOLDER OR THE CODE WILL NOT WORK**

3. Use source activate to activate your qiime2 environment. 
  source activate qiime2-amplicon-2024.5 {this makes qiime2 functions available for use}

4. The first step is to import your sequences with this function (see in line): 
  qiime tools import  \ {this is the function for importing}
   --type 'SampleData[PairedEndSequencesWithQuality]' \ {tells qiime2 what kind of sequences that you're importing}
     --input-path metadata/manifest.csv \ {this directs qiime2 to a manifest file that has the name of your sample, the filepath the sequencing reads from that sample, and the direction 'forward' or 'reverse' from that sequencing run}
      --output-path output/16S_eDNA_Demux.qza \ {this is the name of the output file that contains all of your sequences. }
       --input-format PairedEndFastqManifestPhred33 {this is the format that your sequences are found in.}
       
  qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path metadata/manifest.csv --output-path output/16S_eDNA_Demux.qza --input-format PairedEndFastqManifestPhred33
       
5. Visualize your sequencing reads. Go to QIIME2View.org to view the output *.qzv file. 
  qiime demux summarize \ {summarizes the data you've just imported}
  --i-data output/16S_eDNA_Demux.qza \ {this is the input that you just imported. }
  --o-visualization output/16S_eDNA_Demux.qzv {output is a visualizable file}
  
  qiime demux summarize --i-data output/16S_eDNA_Demux.qza --o-visualization output/16S_eDNA_Demux.qzv
  
  NOTE: Here, you can make sure that the quality of each of your base calls make sense. On the y axis of the interactive quality plot, you'll see 0-40 with 40 meaning 1 error in 10,000 basecalls while 10 means 1 error in every 10 base calls.  If the quality drops below 30 prior to the 200th base, tell Weihua Wang at UIC to see if she'll re-run the samples. 
  
6. Next, primers will be trimmed from these sequences.  As primers represent artificial DNA, we have to remove it prior to downstream sequence variant calling. 

  qiime cutadapt trim-paired \
  --i-demultiplexed-sequences output/16S_eDNA_Demux.qza \
  --p-front-f ACGAGAAGACCCYRTGGARCTT\ {forward primer sequence in forward read}
  --p-front-r ATCCAACATCGAGGTCGTAA\ {reverse primer sequence in reverse read}
  --p-error-rate 0.3\ {higher error rate to allow for primer changes}
  --p-match-adapter-wildcards True \ {this and below allow for non-ACTG bases}
  --p-match-read-wildcards True \
  --p-discard-untrimmed True \ {throws out any read that doesn't have each primer}
  --p-times 10 \ {this is to remove concatenated primers}
  --o-trimmed-sequences output/lead_primertrimmed_16S_eDNA_Demux.qza
  
  qiime cutadapt trim-paired \
  --i-demultiplexed-sequences output/16S_eDNA_Demux.qza \
  --p-front-f ACGAGAAGACCCYRTGGARCTT \
  --p-front-r ATCCAACATCGAGGTCGTAA \
  --p-error-rate 0.3 \
  --p-match-adapter-wildcards True \
  --p-match-read-wildcards True \
  --p-discard-untrimmed True \
  --p-times 10 \
  --o-trimmed-sequences output/lead_primertrimmed_16S_eDNA_Demux.qza
  
  qiime cutadapt trim-paired \
  --i-demultiplexed-sequences output/lead_primertrimmed_16S_eDNA_Demux.qza \
  --p-adapter-f TTACGACCTCGATGTTGGAT\ {reverse complement of reverse primer sequence in forward read}
  --p-adapter-r AAGYTCCAYRGGGTCTTCTCGT\ {reverse complement of forward primer sequence in reverse read}
  --p-error-rate 0.3\ {higher error rate to allow for primer changes}
  --p-match-adapter-wildcards True \ {this and below allow for non-ACTG bases}
  --p-match-read-wildcards True \
  --p-discard-untrimmed False \ 
  --p-times 10 \ {this is to remove concatenated primers}
  --o-trimmed-sequences output/primertrimmed_16S_eDNA_Demux.qza 
  
  qiime cutadapt trim-paired \
  --i-demultiplexed-sequences output/lead_primertrimmed_16S_eDNA_Demux.qza \
  --p-adapter-f TTACGACCTCGATGTTGGAT \
  --p-adapter-r AAGYTCCAYRGGGTCTTCTCGT \
  --p-error-rate 0.3 \
  --p-match-adapter-wildcards True \
  --p-match-read-wildcards True \
  --p-discard-untrimmed False \
  --p-times 10 \
  --o-trimmed-sequences output/primertrimmed_16S_eDNA_Demux.qza 

NOTE: You'll see that we actually performed two primer trimming steps: 1) to remove the forward and reverse primer sequences from the 3' end of the reads, tossing out any thing that doesns't have both of these; and 2) to remove the reverse complements of the forward and reverse primer sequences at the 5' end of the reads, but retaining every read.  This was done to make sure that the primer sequences were removed, but that we did not throw out sequencing reads that might be too long to have both the primer and reverse complement of the primer within the same read.  This is a newly determined step and will not be seen in previous iterations. 

7. Turn this sequence file into a .qzv for visualization on QIIME2.  
  
  qiime demux summarize \
    --i-data output/primertrimmed_16S_eDNA_Demux.qza \
    --o-visualization output/primertrimmed_16S_eDNA_Demux.qzv 
  
  NOTE: Here, you can make sure that the quality of each of your base calls make sense. On the y axis of the interactive quality plot, you'll see 0-40 with 40 meaning 1 error in 10,000 basecalls while 10 means 1 error in every 10 base calls.  
  **MAKE NOTE OF THE FIRST BP WHERE QUALITY DROPS BELOW 30 ON THE FORWARD AND REVERSE READS AS THESE WILL BE USED IN THE FUNCTION BELOW FOR --p-trunc-len-f (forward read) and --p-trunc-len-r (reverse read). 
  
8. Use DADA2 to correct sequencing errors and identify unique sequence variants (ASVs or ESVs depending on who you talk to.). These will be different than OTUs, and these nuances can be found in Hoenig et al. 2021 Review paper in Ornithology as well as others. 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs output/primertrimmed_16S_eDNA_Demux.qza \
  --p-trunc-len-f 242 \ {this uses the quality information above for forward read}
  --p-trunc-len-r 243 \ {this uses the quality information above for forward read
  --p-n-threads 0 \ {this uses all but one core in your computer to speed up the process}
  --output-dir output/denoised_16S_eDNA
  
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs output/primertrimmed_16S_eDNA_Demux.qza \
  --p-trunc-len-f 221 \
  --p-trunc-len-r 207 \
  --p-n-threads 0 \
  --output-dir output/denoised_16S_eDNA

9. Now, we're going turn the .qza files (that can only be used by qiime2) into file types that we can actually work with. 

  qiime tools extract \
  --input-path output/denoised_16S_eDNA/table.qza \
  --output-path output/denoised_16S_eDNA/sample_table

  qiime tools extract \
  --input-path output/denoised_16S_eDNA/representative_sequences.qza \
  --output-path output/denoised_16S_eDNA/representative_sequences

  NOTE: This will make folders that give us the sample table (i.e., which samples have which sequence variants) and the representative sequence variants themselves.  These sequence variants can be uploaded to BLASTN or they can be identified taxonomically using the code below.  However, you will need to make sure you have a 'taxonomic-classification' folder that has the *_seqs.qza and *_taxa.qza files. 
  
10. # Classify Sequences
qiime feature-classifier classify-consensus-blast\
  --i-query output/denoised_16S_eDNA/representative_sequences.qza\
  --i-reference-reads taxonomic-classification/Vertebrata16S_derep1_seqs.qza\
  --i-reference-taxonomy taxonomic-classification/Vertebrata16S_derep1_taxa.qza\
  --p-maxaccepts 50\ {This was a subjective decision and can be modified to allow for more hits to be used in determining the taxonomic classification}
  --p-perc-identity 0.99\ {This was a subjective decision that can be altered to determine how similar a reference sequence must be to be considered for taxonomic classification.}
  --p-query-cov 0.95\ {This was a subjective decision that can be altered to allow more or less of the query sequence to be used in taxonomic classification}
  --output-dir output/denoised_16S_eDNA/classified_taxonomy
  
qiime feature-classifier classify-consensus-blast\
  --i-query output/denoised_16S_eDNA/representative_sequences.qza\
  --i-reference-reads taxonomic-classification/Vertebrata16S_derep1_seqs.qza\
  --i-reference-taxonomy taxonomic-classification/Vertebrata16S_derep1_taxa.qza\
  --p-maxaccepts 50\
  --p-perc-identity 0.99\
  --p-query-cov 0.95\
  --output-dir output/denoised_16S_eDNA/classified_taxonomy

11. Finally we will extract the classified taxonomy which is in qza form into a folder that we can work with.  

12. Proceed to 'r-code/01_making_relevant_files.R' to process the resulting data into something that can be used in data analysis.  

