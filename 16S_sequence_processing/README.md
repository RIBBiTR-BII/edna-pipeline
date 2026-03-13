# 16S Sequence Processing Pipeline

A 16S eDNA amplicon sequence processing workflow
*Created by: Brandon Hoenig & [Cob Staines](https://github.com/cob-staines/)*

## Description

This folder contains Bash (for UNIX-like operating systems) and R scripts to process and analyzing amphibian eDNA results in the [Richards-Zawacki Lab](https://www.rzlab.pitt.edu/), optimized for characterization of amphibian communities from 16S (and some 12S) ribosomal RNA in water samples across the Americas.

This workflow makes assumptions and choices which substantially influence the outputs, some of which are set up as configuration parameters and some of which are scripted. We encourage you to read this documentation, get to know the scripts, and critically observe the workflow parameters and assumptions to ensure that it is appropriate for your application. Suggestions and contributions are welcome!

## Getting Started

### Dependencies

1. Bash shell environment; compatible with Unix-like systems (Linux, macOS) and Windows with Bash support (e.g., WSL, Git Bash).

2. Make sure you have the following installed on your machine:
  - [QIIME2](https://library.qiime2.org/quickstart/amplicon) -- This should be installed using Conda for this workflow (not Docker).
  - You may want to have worked your way through the [Moving Pictures Tutorial](https://docs.qiime2.org/2024.5/tutorials/moving-pictures/) to learn more about this process
  
3. You'll also need to have [R](https://www.r-project.org/) and the following R libraries on your machine:
  -qiime2R
  -tidyverse
  -yaml
  -janitor
  -phylotools
  -readxl
  -openxlsx
  -iucnredlist
  -rgbif
  
Now to the step-by-step guide!
  
### Preparation

1: Clone this repository to a local directory

2: Within the `16S_sequence_processing/runs/` subfolder, create a new run folder with a unique and informative name (e.g. `project_yyyy-mm-dd`). Within your run folder create a folder named `sequences` and copy all Illumina sequence files (`.fastq` or `.fastq.gz`) of interest here (these can be nested in subsequent folders without causing issues).

3: Locate your taxonomic classifier.
  - To build a new NCBI 16S taxonomic classifier, navigate to the `16S_sequence_processing/` folder in the command line and run `bash update-vertebrata16s-classifier.sh`.
  - To use an existing classifier, copy the classifier into a dedicated folder in `16S_sequence_processing/`.

4: Open the `16S_sequence_processing/config.yml` configuration file and adjust to your run specifications following the comments. Save the config file when finished.


### Initiate Sequence Processing
  
To run the workflow, navigate to `16S_sequence_processing/` in the command line, and run `bash run-seq-analysis.sh`. When the workflow is finished, navigate to `your-run-folder/outputs/` to engage with processed outputs for use in subsequent analysis.

## Tools and Tips

 - When you initiate this workflow, the `16S_sequence_processing/config.yml` configuration file is copied to `your-run-dir/output/metadata/config.yml` where it is referenced by the workflow scripts. You can therefore subsequently revise the configuration file and run additional workflows simultaneously, with each run referencing its own local configuration parameters. You the also have a copy of the configuration file to reference and rerun the analysis with ease.
 - When the workflow terminates, a run metadata file `your-run-dir/output/metadata/run_metadata.yml` is generated which contains all of the configuration parameters, as well as some run statistics and diagnostics.

