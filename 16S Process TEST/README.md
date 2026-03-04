# eDNA Sequence Analaysis Pipeline

*Created by: Brandon Hoenig & [Cob Staines](cobstainesconsulting@gmail.com)*
A 16S eDNA amplicon sequence processing workflow

## Description

This repository contains bash and R scripts (for UNIX-like operating systems) to process and analyzing eDNA results in the [Richards-Zawacki Lab](https://www.rzlab.pitt.edu/), optimized for characterization of amphibian communities from 16S ribosomal RNA in water samples across the Americas.

This workflow makes assumptions and choices which substantially influence the outputs, some of which are set up as configuration parameters and some of which are scripted. We encourage you to read this documentation, get to know the scripts, and critically observe the workflow parameters and assumptions to ensure that it is appropriate for your application. Suggestions and contributions are welcome!

## Getting Started

### Dependencies

1. Bash shell environment; compatible with Unix-like systems (Linux, macOS) and Windows with Bash support (e.g., WSL, Git Bash).

2. Make sure you have the following installed on your machine:
  - [QIIME2](https://library.qiime2.org/quickstart/amplicon) -- This should be installed unish Conda for this workflow (not Docker).
  - You may want to have worked your way through the [Moving Pictures Tutorial](https://docs.qiime2.org/2024.5/tutorials/moving-pictures/) to learn more about this process
  
3. You'll also need to have R installed and have these libraries on your machine:
  qiime2R
  tidyverse
  yaml
  janitor
  phylotools
  readxl
  openxlsx
  iucnredlist
  rgbif
  
Now to the step-by-step guide!
  
**Preparation** 

1: Clone this repository to a local directory

2: Within the `runs` subfolder (found in the root folder of this repository), create a run folder with a unique name. Create a folder named `sequences` in your run folder, and copy all Illumina sequence files (FASTQ) of interest here (these can be nested in subsequent folders without issues).

3: Identify your taxonomic classifier. To build a new taxonomic classifier, open the repository root folder in terminal and run `bash update-vertebrata16s-classifier.sh`. Alternatively, copy an existing classifier to an identifiable directory in the root folder.

4: Open the config.yml file and adjust to your needs following the comments. Save the config file when finished.


**Sequence Processing** 
  
1. Open the root folder in terminal, and run `bash run-seq-analysis.sh

2. Navigate to `your-run-folder/outputs/` to engage with processed outputs for use in subsequent analysis.

