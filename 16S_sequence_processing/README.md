# Amphibian 16S Sequence Processing Pipeline

An amphibian 16S eDNA amplicon sequence processing workflow

*Created by: Brandon Hoenig & [Cob Staines](https://github.com/cob-staines/)*

## Description

This folder contains Bash (for UNIX-like operating systems) and R scripts to process and analyzing amphibian eDNA results in the [Richards-Zawacki Lab](https://www.rzlab.pitt.edu/), optimized for characterization of amphibian communities from 16S (and some 12S) ribosomal RNA in water samples across the Americas.

This workflow makes assumptions and choices which substantially influence the outputs, some of which are set up as configuration parameters and some of which are scripted. We encourage you to read this documentation, get to know the scripts, and critically observe the workflow parameters and assumptions to ensure that it is appropriate for your application. Suggestions and contributions are welcome!

## Getting Started

### Setup

1. Install Docker -- Choose one of the options below:
  
  - [Docker Desktop (Windows/Mac/Linux)](https://docs.docker.com/get-started/introduction/get-docker-desktop/): Docker Desktop establishes a dedicated Linux virtual machine on which to run Docker Engine.
  - [Docker Engine (Linux)](https://docs.docker.com/engine/install/): If you already have a Linux machine or virtual machine (e.g. WSL), you can install Docker Engine there directly without having to establish a new Linux virtual machine.

  Confirm your Docker installation in the command line (terminal/PowerShell/Command Prompt.)  with: `docker --version`
  Should return: `Docker version xx.x.x, build xxxxxxx`

  *Note: If you installed Docker Engine within a Linux virtual machine, you will need to run all subsequent steps within this virtual machine.* 

2: [Download](https://github.com/RIBBiTR-BII/edna-pipeline/archive/refs/heads/main.zip) or clone (`https://github.com/RIBBiTR-BII/edna-pipeline.git`) this repository to a local directory which is accessible to your Bash shell.

3: Navigate to the `16S_sequence_processing` subfolder in the command line (e.g. `cd your/path/to/16S_sequence_processing`) 

4: Build the initial temporary docker environment with: `docker build -t edna-pipeline .` This will create a local environment for the workflow following the `dDckerfile` instructions. This will require an active internet connection. It will also take some time the first time you build, as the image must install QIIME, R, and a few other dependencies. A wired internet connection will speed this up.

### Run Sequence Processing

1. Navigate to the `16S_sequence_processing` subfolder in the command line, and build the Docker environment with `docker build -t edna-pipeline .` *This only needs to be done once at the start of a new session. If you are preceding from the Setup steps above, skip to Step 2..*

2. Locate sequences for your run, electing one of the following
  
  - To test the pipeline, an existing run folder `16S_sequence_processing/runs/test_run_01` has been included which contains test sequences. Move on to Step 3.
  - To run the pipeline on your own sequences, create a run folder (e.g. `16S_sequence_processing/runs/run-name_yyy-mm-dd`). Create a `sequences` subfolder (e.g. `16S_sequence_processing/runs/run-name_yyy-mm-dd/sequences`) and copy your .fastq or .fastq.gz amplicon sequence files here. Sequence files can be nested in subfolders without issue.

3. Build a local 16S classifier *(Optional)*: In the command line, run the following command:
  
  - **Windows PowerShell:** `docker run -it --rm --user $(id -u):$(id -g) -v ${PWD}:/data edna-pipeline bash build-amphibia16s-classifier.sh`
  - **Windows Command Prompt:** `docker run -it --rm --user $(id -u):$(id -g) -v %cd%:/data edna-pipeline bash build-amphibia16s-classifier.sh`
  - **macOS/Linux:** `docker run -it --rm --user $(id -u):$(id -g) -v $(pwd):/data edna-pipeline bash build-amphibia16s-classifier.sh`
  
  You can alternatively copy your own classifier to the `16S_sequence_processing` folder, or skip using a local classifier.

4. Open `config.yml` in a text browser, and confirm the configuration settings for your run. Save this file after editing. By default, `config.yml` is set up for the test run `test_run_01` without a local classifier. To run your own sequences, specify the run directory containing your sequences.

5. Run sequence processing:
  
  - **Windows PowerShell:** `docker run -it --rm --user $(id -u):$(id -g) -v ${PWD}:/data edna-pipeline bash run-seq-processing.sh`
  - **Windows Command Prompt:** `docker run -it --rm --user $(id -u):$(id -g) -v %cd%:/data edna-pipeline bash run-seq-processing.sh`
  - **macOS/Linux:** `docker run -it --rm --user $(id -u):$(id -g) -v $(pwd):/data edna-pipeline bash run-seq-processing.sh`

## Tools and Tips

 - When you initiate this workflow, the configuration file (`16S_sequence_processing/config.yml`) is copied to `[YOUR RUN DIR]/output/metadata/config.yml` where it is referenced by the workflow scripts. You can therefore subsequently revise the configuration file and run additional workflows simultaneously, with each run referencing its own local configuration parameters. You then also have a copy of the configuration file to reference and rerun the analysis in the future.
 - When the workflow terminates, a run metadata file `your-run-dir/output/metadata/run_metadata.yml` is generated which contains all of the configuration parameters, as well as some run statistics and diagnostics.

