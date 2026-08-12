# Analysis of Prcoessed 16S Sequences Results from the 

A workflow for analyzing results from the [Amphibian 16S Sequence Processing Pipeline](https://github.com/RIBBiTR-BII/edna-pipeline/tree/main/16S_sequence_processing)

*Created by: Brandon Hoenig & [Cob Staines](https://github.com/cob-staines/)*

## Description
At this point you have successfully processed some sequences through the [Amphibian 16S Sequence Processing Pipeline](https://github.com/RIBBiTR-BII/edna-pipeline/tree/main/16S_sequence_processing). Great work! Now what? This folder contains a series of R scripts to begin to analyze your results. These steps include taxonomic classification of ASVs, aligning samples with field survey data, and controlling for contamination. Follow the steps below to begin this analysis.

## Setup
1. **Create an Entrez API key:** This will allow you to look up NCBI taxa through the `taxize` package without additional steps (used in `analysis/general/r/04_web_blast_json_parse.Rmd`):

    - [Register with NCBI](https://account.ncbi.nlm.nih.gov/signup/) (if you have not already)
    - Log in and generate an Entrez API key following [this guidance](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317). Copy your API key.
    - Save this key to your local `.Renviron`. In RStudio, run the following lines:
    
      ```{r}
      install.packages("usethis")  # if not already installed
      usethis::edit_r_environ()
      ```
  
      In the .Renviron document that opens, save your copied API key as: 
  
          `ENTREZ_KEY = "your-key-here".`
          
      Save and close the .Renviron document. Then in the RStudio menu go to `Session` -> `Restart R`. You can test that your key is saved and accessible by running:
 
      ```{r}
      Sys.getenv("ENTREZ_KEY")
      ```
      
      This should display your Entrez API key (an empty string `""` means this was unsuccessful).

    

2. **Set up GBIF API access:** This will allow you to look up GBIF occurrences of taxa of interest without additional steps (used in `analysis/general/r/05_query_taxonomy_geography.Rmd`):

    - [Register with GBIF](https://www.gbif.org/) -> Login *(top right)* -> Register (If you have not already)
    - Follow [this tutorial](https://docs.ropensci.org/rgbif/articles/gbif_credentials.html) to save your GBIF login where it can be accessed by the package `rgbif`.


3. **Establish a connection to RIBBiTR database:** This will allow you to link eDNA results with collection metadata (used in `analysis/general/r/07_join_sample_collection.Rmd`):

    - Follow the [RIBBiTR DB connection tutorial](https://ribbitr-bii.github.io/ribbitr-data-access/tutorial_series/01_connection_setup.html) to connect in RStudio.

## Analysis
These scripts have not been automated, as they contain more decisions for users to consider as they conduct the analysis. To begin, navigate to the `analysis/general/r/` folder. Each script here corresponds to the steps below. Open each scrip in RStudio, look at the header notes, set the configuration ("config") sections to meet your needs and run the script:

4. **Web Blast & Parse** *(`04_web_blast_json_parse.Rmd`)*: Follow script instructions to upload the representative sequences to [NCBI's Web Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) service, and download the .json query results. Run the script to parse the .json outputs.

5. **GBIF Query** *(05_query_taxonomy_geography.Rmd)*: This script searches for occurrences of reference taxonomies in the study system of interest, to prioritize classification of local species. This pulls in hits from any of the following sources: BLAST, Vsearch, or Web BLAST.

6. **Classify ASVs** *(06_classify_asv.Rmd)*: This script uses a hierarchy of classification methods to assign taxonomic hits to each ASV, pulling from the GBIF query (if conducted) and incorporating hits from any of the following sources: : BLAST, Vsearch, or Web BLAST.

7. **Map Samples** *(07_sample_map.Rmd)*: This script Maps Illumina samples to RIBBiTR sample ids, to support alignment of results with collection metadata downstream.

8. **Control for Contamination** *(08_sample_controls.Rmd)*: This script corrects ASV read numbers for contamination by removing up to 100x the read numbers of any ASVs found in associated controls.

9. **Export Results** *(09_export_results.Rmd)*: This script combines results from steps 6, 7, and 8, and as well as sample metadata from thge RIBBiTR database, to create two cohesive outputs: 1) for ASVs (reads, classifications, etc.) and 2) field samples (collection site, date, filter method, etc.)
