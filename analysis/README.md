# Analysis of Prcoessed 16S Sequences Results from the 

A workflow for analyzing results from the [Amphibian 16S Sequence Processing Pipeline](https://github.com/RIBBiTR-BII/edna-pipeline/tree/main/16S_sequence_processing)

*Created by: Brandon Hoenig & [Cob Staines](https://github.com/cob-staines/)*

## Description
At this point you have successfully processed some sequences through the [Amphibian 16S Sequence Processing Pipeline](https://github.com/RIBBiTR-BII/edna-pipeline/tree/main/16S_sequence_processing). Great work! Now what? This folder contains a series of R scripts to begin to analyze your results. These steps include taxonomic classification of ASVs, aligning samples with field survey data, and controlling for contamination. Follow the steps below to begin this analysis.

## Setup
1. **Create an Entrez API key:** This will allow you to look up NCBI taxa through the `taxize` package without additional steps (used in `analysis/general/r/04_web_blast_json_parse.Rmd`):

    - [Register with NCBI](https://account.ncbi.nlm.nih.gov/signup/) (if you have not already)
    - Log in and generate an Entrez API key following [this guidance](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317). Copy your API key.
    - In your RStudio Console, run the following lines to save the API key to your local `.Renviron`:
    
      ```{r}
      install.packages("usethis")
      usethis::edit_r_environ()
      ```
  
      In the .Renviron document that opens, save your copied API key as: 
  
      `ENTREZ_KEY = "your-key-here"`
          
      Save and close the .Renviron document. Then in the RStudio menu go to `Session` -> `Restart R`. You can test that your key is saved and accessible by running:
 
      ```{r}
      Sys.getenv("ENTREZ_KEY")
      ```
      
      This should display your Entrez API key (an empty string `""` means the key is not found).

2. **Set up GBIF API access (optional):** This will allow you to look up GBIF occurrences of taxa (used optionally in `analysis/general/r/05_query_taxonomy_geography.Rmd`):

    - [Register with GBIF](https://www.gbif.org/) -> Login *(top right)* -> Register (If you have not already)
    - Follow [this tutorial](https://docs.ropensci.org/rgbif/articles/gbif_credentials.html) to save your GBIF login where it can be accessed by the package `rgbif`.

3. **Establish a connection to RIBBiTR database:** This will allow you to link eDNA results with field collection metadata (used in `analysis/general/r/07_join_sample_collection.Rmd`):

    - Follow the [RIBBiTR DB connection tutorial](https://ribbitr-bii.github.io/ribbitr-data-access/tutorial_series/01_connection_setup.html) to connect in RStudio.

4. **Create an RStudio project (optional):** Open RStudio, select `File -> New Project -> Existing Directory -> Browse` and browse to your local directory of this `edna-pipeline` repository. Thel select `Create Project`. This is not required, but will make it easier to navigate between the various analysis scripts.

## Analysis
The numbered (4 - 9) analysis steps below correspond to numbered .Rmd scripts which should be run in RStudio in succession. They have not been automated, as each script contains decisions for users to consider as the analysis progresses. To begin, navigate to the `analysis/general/r/` folder. Open each script in RStudio, review the header notes, set the parameters in the configuration ("config") sections at the top to meet your needs, and run each script.

4. **Web Blast & Parse** *(`04_web_blast_json_parse.Rmd`)*: Follow script instructions below to upload the representative sequences to [NCBI's Web Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) service, and download the query results. This script parses the .json outputs from the Web BLAST query.
    - a. Upload the ASV representative sequences .fasta file to NCBI's Web BLAST: Nucleotide BLAST service
      - Visit https://blast.ncbi.nlm.nih.gov/Blast.cgi, click on Nucleotide BLAST
      - In the `Enter Query Sequence` panel, beside `Or, upload file`, click `Browse` and navigate to the ASV representative sequences .fasta file at: `[your-run-directory]/analysis/s06_denoised_16S_eDNA/representative sequences/.../data/dna-sequences.fasta`
      - Add a descriptive `Job Title`
      - Under `Program Selection: Optimize for`, select `More dissimilar sequences (discontiguous megablast)` (ideal for eDNA)
      - Under `Algorithm parameters`
        - adjust the max number of hits as desired (10 - 100 is likely fine)
        - adjust the `Expect threshold` (0.03 is recommended)
      - Click `BLAST` and wait for the query to finish
    - b. In the main Web BLAST results panel, to the right of `RID`, click `Download All` and select `Single-file JSON`. Save the JSON report file to `[your-run-directory]/outout/`
    - c. Once you have the results file, adjust the parameters in the `Config` section to match your needs. You can then run this script to parse and structure the results for downstream analysis.

5. **GBIF Query** *(05_query_taxonomy_geography.Rmd)*: This script searches for occurrences of reference taxonomies in the study system of interest, to prioritize classification of local species and provide context for interpretation. This pulls in hits from any of the following sources: BLAST, Vsearch, or Web BLAST.
  - This script requires API keys for GBIF in .Renviron (see `Setup` above).
  - This step is optional. If you want to skip this step, proceed to the next script (`06_classify_asv.Rmd`) and set config the parameter `gbif_query` to `FALSE`.
  
6. **Classify ASVs** *(06_classify_asv.Rmd)*: This script uses a hierarchy of classification methods to assign taxonomic hits to each ASV, optionally pulling from the GBIF query and incorporating hits from any of the following sources: : BLAST, Vsearch, or Web BLAST.
  - A likely taxonomy is assigned to each ASV following the hierarchy `accept_method`s if the given `accept_method` criteria are met. All assigned taxonomy from all methods, along with all hits, are exported to the specified `hybrid_classification_out` path.
  - The single "best" classifications (i.e. `accept_method` with greatest priority in hierarchy) for each ASV are exported to the specified `classification_out` path.


7. **Map Samples** *(07_sample_map.Rmd)*: This script Maps Illumina samples to RIBBiTR sample ids, to support alignment of results with collection metadata downstream.\
  - This requires a connection to the RIBBiTR database (see `Setup` above).

8. **Control for Contamination** *(08_sample_controls.Rmd)*: This script calculates controlled ASV reads for each sample by subtracting any read counts found in controls (multiplied by a scaling factor `asv_control_th_factor`), with minimum controlled reads of 0.
  - Lab positive and negative controls are applied to all samples globally, while field negative controls are applied to corresponding field samples only.

9. **Export Results** *(09_export_results.Rmd)*: This script combines results from steps 6, 7, and 8, and as well as sample metadata from the RIBBiTR database, to create two cohesive outputs:
  a. for ASVs (reads, classifications, etc.)
  b. field samples (collection site, date, filter method, etc.)
  
## After Preliminary Analysis

You are now ready to move into your own analysis to address the research questions you have at hand!

**Note on archiving run directories:** Once you are confident that you are done with the processing and analysis workflows for a given dataset, you may consider deleting the `sequences` and `analysis` folders from your run directory. The remaining `outputs` folder should have all the outputs from theis pipeline, and metadata on the sequence files used and the configuration settings if needed for future reference. This will help you reclaim some space on your machine (just make sure that the sequence files are archived somewhere else!).
