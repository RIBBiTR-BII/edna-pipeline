# Amphibian 16S Sequence Processing Pipeline

An amphibian 16S eDNA amplicon sequence processing workflow

*Created by: Brandon Hoenig & [Cob Staines](https://github.com/cob-staines/)*

## Description

This folder contains Bash (for UNIX-like operating systems) and R scripts to process and analyzing amphibian eDNA results in the [Richards-Zawacki Lab](https://www.rzlab.pitt.edu/), optimized for characterization of amphibian communities from 16S (and some 12S) ribosomal RNA in water samples across the Americas.

This workflow makes assumptions and choices which substantially influence the outputs, some of which are set up as configuration parameters and some of which are scripted. We encourage you to read this documentation, get to know the scripts, and critically observe the workflow parameters and assumptions to ensure that it is appropriate for your application. Suggestions and contributions are welcome!

## Getting Started

### Setup

1. Install Podman -- Choose one of the options below:

   - **Windows:** Install [Podman Desktop for Windows](https://podman-desktop.io/). This will also set up Podman Machine, a lightweight Linux virtual machine (VM) that Podman uses on Windows. During first launch, follow the prompt to initialise Podman Machine.
   - **macOS:** Install [Podman Desktop for Mac](https://podman-desktop.io/). After installation, initialise Podman Machine from the Podman Desktop interface or run `podman machine init && podman machine start` in the terminal.
   - **Linux:** Install Podman Engine directly — no VM required:
     ```bash
     sudo apt-get update && sudo apt-get install -y podman   # Debian/Ubuntu
     ```

   Confirm your Podman installation in the command line (terminal/PowerShell/Command Prompt) with: `podman --version`  
   Should return: `podman version x.x.x`

   > **Docker alternative:** This workflow is also fully compatible with [Docker](https://docs.docker.com/get-started/introduction/get-docker-desktop/). Simply substitute `podman` for `docker` in all commands below.

   > **Note for WSL2 users:** If you are running a Linux virtual machine on Windows via WSL2, you can install Podman Engine inside the VM directly and run all commands from within the VM terminal, without needing Podman Desktop.

2. [Download](https://github.com/RIBBiTR-BII/edna-pipeline/archive/refs/heads/main.zip) or clone (`git clone https://github.com/RIBBiTR-BII/edna-pipeline.git`) this repository to a local directory accessible to your command line.

3. Navigate to the `16S_sequence_processing` subfolder in the command line:
   ```bash
   cd your/path/to/16S_sequence_processing
   ```

4. Build the container image:
   ```bash
   podman build -t edna-pipeline .
   ```
   This creates a local container image following the `Dockerfile` instructions, installing QIIME2, R, and all required dependencies. An active internet connection is required. **This only needs to be done once** — or until the Dockerfile changes. A wired connection will speed this up considerably as the base image is several gigabytes.

---

### Test Run

1. Navigate to the `16S_sequence_processing` subfolder in the command line. If you have not yet built the container image, run `podman build -t edna-pipeline .` first (see Setup Step 4 above).

2. Run sequence processing:

   - **Windows PowerShell:**
     ```powershell
     podman run -it --rm -v ${PWD}:/data:z edna-pipeline bash run-seq-processing.sh
     ```
   - **Windows Command Prompt:**
     ```cmd
     podman run -it --rm -v %cd%:/data:z edna-pipeline bash run-seq-processing.sh
     ```
   - **macOS/Linux:**
     ```bash
     podman run -it --rm -v $(pwd):/data:z edna-pipeline bash run-seq-processing.sh
     ```

---

### Sequence Processing Run

1. Navigate to the `16S_sequence_processing` subfolder in the command line. If you have not yet built the container image, run `podman build -t edna-pipeline .` first (see Setup Step 4 above).

2. Build a local Amphibia 16S classifier: Run the following command:

   - **Windows PowerShell:**
     ```powershell
     podman run -it --rm -v ${PWD}:/data:z edna-pipeline bash build-amphibia16s-classifier.sh
     ```
   - **Windows Command Prompt:**
     ```cmd
     podman run -it --rm -v %cd%:/data:z edna-pipeline bash build-amphibia16s-classifier.sh
     ```
   - **macOS/Linux:**
     ```bash
     podman run -it --rm -v $(pwd):/data:z edna-pipeline bash build-amphibia16s-classifier.sh
     ```

   You can alternatively copy your own classifier to the `16S_sequence_processing` folder, or skip this step to run without a local classifier.

3. Locate sequences for your run, selecting one of the following:

   - **Test the pipeline:** An existing run folder `16S_sequence_processing/runs/test_run_01` has been included with test sequences. Proceed to Step 3.
   - **Your own sequences:** Create a run folder (e.g. `16S_sequence_processing/runs/run-name_yyyy-mm-dd`), create a `sequences` subfolder inside it, and copy your `.fastq` or `.fastq.gz` paired amplicon sequence files there. Sequence files can be nested in subfolders without issue.


4. Open `config.yml` in a text editor and confirm the configuration settings for your run. Save the file after editing. By default, `config.yml` is set up for the test run `test_run_01` without a local classifier. To run your own sequences, update the run directory and any other relevant parameters.

5. Run sequence processing:

   - **Windows PowerShell:**
     ```powershell
     podman run -it --rm -v ${PWD}:/data:z edna-pipeline bash run-seq-processing.sh
     ```
   - **Windows Command Prompt:**
     ```cmd
     podman run -it --rm -v %cd%:/data:z edna-pipeline bash run-seq-processing.sh
     ```
   - **macOS/Linux:**
     ```bash
     podman run -it --rm -v $(pwd):/data:z edna-pipeline bash run-seq-processing.sh
     ```

---

## Tools and Tips

- When you initiate this workflow, the configuration file (`16S_sequence_processing/config.yml`) is copied to `[YOUR RUN DIR]/output/metadata/config.yml` where it is referenced by the workflow scripts. You can therefore subsequently revise the configuration file and run additional workflows simultaneously, with each run referencing its own local configuration parameters. You also have a copy of the configuration file to reference and rerun the analysis in the future.

- When the workflow terminates, a run metadata file `your-run-dir/output/metadata/run_metadata.yml` is generated which contains all configuration parameters, as well as run statistics and diagnostics.

- **File permissions:** Podman runs containers as your own user by default (rootless), so output files and folders created by the pipeline will be owned by you and fully accessible after the run. If you are using Docker instead and encounter permission issues on output folders, add `--user $(id -u):$(id -g)` to your `docker run` commands, or reclaim ownership of existing folders with `sudo chown -R $USER:$USER runs/`.

- **The `:z` volume flag** on the `-v` mount option sets the correct SELinux label on Linux systems, allowing the container to read and write your project files. It is safely ignored on macOS and Windows.

- **Debugging interactively:** To open a shell inside the container for manual testing or troubleshooting, run:
  ```bash
  podman run -it --rm -v $(pwd):/data:z edna-pipeline bash   # macOS/Linux
  ```
  Your project files will be visible at `/data`. Type `exit` to leave.