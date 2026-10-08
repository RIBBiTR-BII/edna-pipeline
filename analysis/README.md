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

3. **Establish a connection to RIBBiTR database:** This will allow you to link eDNA results with field collection metadata (used in `analysis/general/r/03_sample_map.Rmd`):

    - Follow the [RIBBiTR DB connection tutorial](https://ribbitr-bii.github.io/ribbitr-data-access/tutorial_series/01_connection_setup.html) to connect in RStudio.

4. **Create an RStudio project (optional):** Open RStudio, select `File -> New Project -> Existing Directory -> Browse` and browse to your local directory of this `edna-pipeline` repository. Thel select `Create Project`. This is not required, but will make it easier to navigate between the various analysis scripts.

## Analysis
The numbered (3 - 9) analysis steps below correspond to numbered .Rmd scripts which should be run in RStudio in succession. To begin, navigate to the `analysis/general/r/` folder. (Step 7 is reserved for an upcoming per-ASV home-system / contamination-detection script and doesn't exist yet -- the sequence currently skips from 6 to 8.)

Before running the pipeline, open `00_pipeline_config.yml` and review/update its parameters (run directory, thresholds, etc.) to match your run. All scripts 03-06 and 08-09 read this shared config file, so it only needs to be edited once per run. A run may contain samples from multiple study systems -- which system(s) are present is derived automatically from the sample map (step 3), not configured by hand.

You then have two options for running the scripts:
- **Step through manually:** Open each script in RStudio, review the header notes, and run it chunk by chunk. This is recommended the first time through, as each script contains decisions for users to consider as the analysis progresses.
- **Run end to end:** Once you're comfortable with the decisions each script makes, source `00_run_pipeline.R` to render scripts 03-06 and 08-09 in sequence using the settings in `00_pipeline_config.yml`.

3. **Map Samples** *(`03_sample_map.Rmd`)*: This script maps Illumina samples to RIBBiTR sample ids and assigns each sample a `study_system`, to support alignment of results with collection metadata downstream and let later steps (starting with step 5) know which study system(s) are present in the run.
  - This requires a connection to the RIBBiTR database (see `Setup` above).
  - Runs first, since step 5 depends on its `study_system` column.

4. **Web Blast & Parse** *(`04_web_blast_json_parse.Rmd`)*: Follow script instructions below to upload the representative sequences to [NCBI's Web Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) service, and download the query results. This script parses the .json outputs from the Web BLAST query.
    - a. Upload the ASV representative sequences .fasta file to NCBI's Web BLAST: Nucleotide BLAST service
      - Visit https://blast.ncbi.nlm.nih.gov/Blast.cgi, click on Nucleotide BLAST
      - In the `Enter Query Sequence` panel, beside `Or, upload file`, click `Browse` and navigate to the ASV representative sequences .fasta file at: `[your-run-directory]/analysis/s06_denoised_16S_eDNA/representative sequences/.../data/dna-sequences.fasta`
      - Add a descriptive `Job Title`
      - Under `Program Selection: Optimize for`, select `More dissimilar sequences (discontiguous megablast)` (ideal for eDNA)
      - Under `Algorithm parameters`
        - adjust the max number of hits as desired (50 is likely fine)
        - adjust the `Expect threshold` (0.03 is recommended)
      - Click `BLAST` and wait for the query to finish
    - b. In the main Web BLAST results panel, to the right of `RID`, click `Download All` and select `Single-file JSON`. Save the JSON report file to `[your-run-directory]/outout/`
    - c. Once you have the results file, adjust the parameters in the `Config` section to match your needs. You can then run this script to parse and structure the results for downstream analysis.

5. **GBIF Query** *(05_query_taxonomy_geography.Rmd)*: This script searches for occurrences of reference taxonomies in each study system present in the run (per step 3's sample map), to prioritize classification of local species and provide context for interpretation. This pulls in hits from any of the following sources: BLAST, Vsearch, or Web BLAST.
  - This script requires API keys for GBIF in .Renviron (see `Setup` above).
  - This step is optional. If you want to skip this step, proceed to the next script (`06_classify_asv.Rmd`) and set config the parameter `gbif_query` to `FALSE`.
  - Taxon resolution runs once per run; occurrence counts are queried and checkpointed separately per study system, so the script can be re-run to resume after a partial failure or once more systems' samples are added.

6. **Classify ASVs** *(06_classify_asv.Rmd)*: This script uses a hierarchy of classification methods to assign taxonomic hits to each ASV, optionally pulling from the GBIF query and incorporating hits from any of the following sources: : BLAST, Vsearch, or Web BLAST.
  - A likely taxonomy is assigned to each ASV following the hierarchy `accept_method`s if the given `accept_method` criteria are met. All assigned taxonomy from all methods, along with all hits, are exported to the specified `hybrid_classification_out` path.
  - The single "best" classifications (i.e. `accept_method` with greatest priority in hierarchy) for each ASV are exported to the specified `classification_out` path.

8. **Contamination Filtering** *(08_contamination_filtering.Rmd)*: This script turns the flags from steps 6 and 7, and the read counts in the run's controls, into the clean read count for every sample x ASV. It is the only step that removes reads.
  - Each row keeps the raw count next to the clean count; a removed detection gets a clean count of 0 and a `filter_reason` (`short_sequence`, `contaminant_library`, `anthropogenic`, `positive_control`, or `control_threshold`). Each filter can be switched off under `filters` in `00_pipeline_config.yml`.
  - Control threshold: PCR and extraction negative controls are applied to all samples globally, while field negative controls are applied to corresponding field samples only (raw count minus `asv_control_th_factor` x the control's count, floored at 0).
  - PCR positive-control components are identified by sequence from the positive controls, and removed outside the systems where they are local.

9. **Export Results** *(09_export_results.Rmd)*: This script combines results from steps 3, 6, and 8, and as well as sample metadata from the RIBBiTR database, to create two cohesive outputs:
  a. for ASVs (reads, classifications, etc.)
  b. field samples (collection site, date, filter method, etc.)
  
## Exploratory: Contamination Classification (steps 11-14)
These scripts develop a classification of every detection as signal or contamination, as a possible successor to step 8's filters. They are not run by `00_run_pipeline.R`, and nothing downstream reads their outputs. Their settings are under `network_*`, `known_contaminants`, `edge_classification` and `validation` in `00_pipeline_config.yml`.

11. **Network Analysis** *(11_network_analysis.Rmd)*: builds the sample--ASV network, exports it for Gephi / Cytoscape, and tests ASV and sample communities for association with lab batches vs. ecology.
12. **ASV Fingerprints** *(12_asv_fingerprints.qmd)*: an interactive viewer (one row of panels per ASV: field samples and field negatives per system, then lab controls) for spotting contamination patterns by eye. It shows raw detections only, with no step 8 or model output; clicking an ASV's id or length copies the full id or its sequence. It also holds the blind labelling mode for the validation sample. Render it with `quarto render` (or the Render button in RStudio).
13. **Edge Classification** *(13_edge_classification.Rmd)*: computes the features of every detection in a field library, fits the model ladder below on the development labels, classifies every detection with each model and baseline, and scores them (only when `validation$report_performance` / `evaluate_test_set` are switched on). It also writes `<run>_13_network_edges_detection.csv`: every detection edge with each model's decision (`control` for control libraries) and probability of signal, to load into Gephi / Cytoscape after script 11's nodes table and colour or filter the network by model.
    - *13_legacy_edge_classification.Rmd*: the legacy iterated log-odds model (v0-v2.1), kept as the `legacy` baseline; run it before script 13 to include that baseline.
14. **Validation Sample** *(14_validation_sample.Rmd)*: draws the seeded, stratified random sample of detections that is labelled by hand and used to score models.

### Validation protocol
*Agreed 2026-10-05, revised 2026-10-07 before any labelling. Follow it for every model change, so that models are compared scientifically rather than tuned to individual cases.*

1. **Reference standard.** A random sample of detections, labelled by hand as `signal`, `contamination` or `uncertain`. These labels are the reference standard. Hand-picked reference cases (`13_edge_labels.csv`) are sanity checks only, and are never used to score or tune a model.
2. **Label definitions.**
    - **Signal**: DNA that was in the environment at the sampled site when the sample was taken.
    - **Contamination**: DNA introduced at any later step -- field handling and gear carryover between sites, extraction, PCR, or sequencing (cross-talk). DNA carried on gear from another site was in the sample tube, but is contamination.
    - The label is about where the DNA came from, not which ASV it is: the same ASV can be signal in one library and contamination in another.
    - DNA can be physically in the sample and still be contamination. Human DNA is the clearest case: from the environment (e.g. sewage, people upstream) it is signal; from the sampler's hands it is contamination. Field negatives receive the same handling as samples, so DNA at levels comparable to the field negatives is consistent with handling. Use `uncertain` when the pattern cannot tell the sources apart.
    - Controls (never sampled, conventions for the models): every detection in a negative control is contamination; in a positive control the spiked DNA is the expected content and everything else is contamination.
3. **Known contaminants.** ASVs that are contamination wherever they appear, by knowledge from outside this data (lab stocks, positive-control constructs), are listed explicitly in `known_contaminants`, each with a reason. The list is hand-reviewed, with no taxon or length rules, so a genuine field taxon is never swept up (human DNA, for example, may be real signal in field water). Known contaminants are clamped as contamination in every model and left out of the validation population: they cannot be validated, so the bar for listing an ASV is high. The 2026-10-07 list (36 ASVs: the Oophaga lab stock and two positive-control constructs with their error variants) was reviewed from a candidate table of ASVs similar to the construct or Oophaga and amphibian ASVs < 120 bp; the marker's natural amplicon is ~225-267 bp for these taxa.
4. **Population.** Every detection (raw count > 0) in a field library, excluding ASVs hard-flagged in step 6 (short sequence, contaminant library) and the known contaminants. The population is fixed by the config, never by a model's output, so it does not move with the model; metrics hold for this population. Control libraries are not sampled.
5. **Strata.** Group (`amphibian` = step 6 class Amphibia, `other`) x study system x abundance band (x = reads / library reads excluding hard flags; bands < 1e-3, 1e-3 to 1e-1, >= 1e-1). Strata are defined from the data alone, never from any model's output, so the same sample is fair to every model.
6. **Allocation.** 20 detections per amphibian stratum and 10 per other stratum, drawn with a fixed seed; a stratum smaller than that is taken whole. Amphibians are oversampled because they are the research target; weighting keeps estimates unbiased for both domains.
7. **Development / test split.** Within each stratum, 2/3 of the draw goes to the **development** set and 1/3 to the **test** set (a fractional detection is assigned at random), before any labelling. Models are fitted and compared on the development set. The test set stays **sealed** (`validation$evaluate_test_set: false`) and is scored only once, for the final comparison of 1-3 candidate models; record that scoring in the change log.
8. **Blind labelling.** Label in script 12's labelling mode, which shows one sampled detection at a time (the viewer shows no step 8 or model output). Use everything an expert would: the fingerprint pattern, taxonomy, GBIF, lab knowledge (e.g. species and constructs worked with in the lab). Do not browse model output for sampled ASVs before labelling them. Add a short note for the reason where useful, and use `uncertain` when the pattern cannot tell. Labels autosave in the browser; export them with "export labels CSV" to the run's `output/validation/` folder as `<run>_validation_labels.csv`, then re-render scripts 12 and 13. GBIF is shown while labelling and is also a model feature, so its contribution to the models may be somewhat overstated; this is accepted, since an expert would use it.
9. **Costs and decisions.** A false signal (contamination called signal) costs **2**, a false contamination (signal called contamination) costs **1**, and leaving a detection unassigned (for manual review) costs **0.25**. The costs are not used to fit models; they turn a model's probability of signal p into a decision, by choosing the decision of least expected cost: contamination if p < 0.25, signal if p > 0.875, unassigned in between. Changing the costs moves the thresholds and the scores, without refitting.
10. **Metrics.** Each labelled detection is weighted by its stratum's population over the number of labelled detections of that stratum in the split. Detections labelled `uncertain` are neither fitted nor scored. **Primary metric: weighted cost per detection**, reported for the **amphibian** domain (primary, the research objective) and **all eDNA** (secondary), alongside coverage, false-signal rate, false-contamination rate and accuracy among assigned detections.
11. **Fitting and model selection.**
    - Models are logistic regressions on the development labels (`signal` vs `contamination`), with each coefficient's sign fixed in advance from its mechanism and light ridge shrinkage (a Normal(0, `prior_sd`) prior on the coefficients of the standardised features); fitted unweighted, scored weighted.
    - Training and scoring are separate: script 13 always fits the models, but shows development-set performance only when `validation$report_performance` is true, so models can be built without seeing how they do.
    - Models are compared by repeated (20 x) 5-fold cross-validation within the development set, with folds **grouped by ASV x study system**, since detections of the same ASV in a system share most features. Baselines and the legacy model are scored on the same folds, so differences are paired.
    - Choose the **simplest rung within one standard error** of the best rung's cross-validated cost. Check it with a drop-one analysis (remove each mechanism from the top model in turn), so a weak early rung cannot hide a useful later one.
    - Final candidates are scored on the test set, with 95% stratified bootstrap intervals for each model's cost and for paired differences from the baselines.
12. **Rules for changing the model.**
    - Every change is motivated by a contamination mechanism, written as a general rule. No case-by-case adjustments.
    - Write the change and its rationale in the change log **before** seeing its effect on the development set.
    - Development-set errors suggest mechanisms to examine; they are not targets to tune away.
    - Prefer fewer terms and parameters: a change that does not clearly improve the cross-validated cost is not kept.
    - Bump `edge_classification$model_version` with every change.
13. **More labels.** If intervals are too wide, raise `validation$wave` and re-run step 14 to draw more detections from the same strata. A sample with labels is never redrawn or overwritten.
14. **Files.** The sample and labels live in the run's `output/validation/` folder (`validation$dir`, relative to `run_dir`), so each run has its own validation set. Run outputs are not tracked by git: back the labels up with the run's outputs, since they are hand work that cannot be regenerated.

### Model ladder
*Planned 2026-10-07, before labelling.* Each rung adds one contamination mechanism. Most features are properties of an ASV in a study system; the edge itself adds its abundance and its library.

**Shared quantities.**
- x: an edge's reads over its library's reads, excluding hard-flagged and known-contaminant reads.
- mu(a, g): the mean x of ASV a over the libraries of group g, zeros included (detection rate x mean abundance when detected). Groups: each system's field samples, and the lab negatives (extraction and PCR negatives).
- **Home share** = mu(a, this system) / (mu(a, this system) + max over other systems of mu(a, other system)): 1 when the ASV is found only here, 0.5 when as much elsewhere, near 0 for a faint copy of another system's taxon.
- **Kit share** = mu(a, this system) / (mu(a, this system) + mu(a, lab negatives)).
- **Human share**: the share of the rest of an edge's library (its own reads left out) from ASVs of genus *Homo*, as an indicator of low target biomass. It does not claim the human reads are contamination.

**Baselines** (fixed rules, no fitting):

| baseline | rule |
|---|---|
| `all_signal` | everything is signal |
| `script_08` | removed by step 8 = contamination |
| `neg_presence` | contamination if the ASV is detected in any of the library's linked negatives (field, extraction, PCR) -- common eDNA practice |
| `neg_subtract` | subtract the ASV's largest read count in the library's linked negatives; contamination if nothing is left -- common eDNA practice |
| `legacy` | the iterated log-odds model, v2.1 |

**Rungs** (each adds to the one above):

| rung | mechanism | feature | sign |
|---|---|---|---|
| B1 | abundance | log10(x) | + |
| M1 | home system: a dense cloud in one system is signal there and cross-talk / lab carryover elsewhere | home share | + |
| M1b | positive-control source | the ASV is a step 8 positive-control component (replaces clamping them) | - |
| M2 | GBIF locality | local (incl. genus only); not local (other system, not in run); anthropogenic and unassessable are neutral | + / - |
| M3 | kitome: as abundant in lab negatives as in field samples | kit share | + |
| M4 | contaminated library: faint detections in low-biomass libraries are suspect, and clean libraries do not vouch for anything | the edge's human share | - |
| M5 | ASV found mostly in dirty samples | the mean human share of the ASV's other edges in the system, minus the mean over all edges of the system (0 for an ASV with no other edge there) | - |

The top model has about 9 fitted parameters. Labels of neighbouring edges are not used as features (no iterative propagation): the neighbour information enters through M4 and M5, from data alone.

**Refinements**, only if the chosen rung leaves errors they would explain, logged first: run-aware home share (other systems sharing a sequencing run only, separating tag-jumping from lab carryover); within-system trace (log10 of x over the ASV's largest x in the system); field negatives (field or gear carryover); detection rate and abundance as separate features; an amphibian intercept.

### Model change log
| version | date | change | rationale | development cost (amphibian / all eDNA) |
|---|---|---|---|---|
| v0 | 2026-10-03 | Initial model: U1 own negatives, U2 cross-talk, U3 replicates, U4 control share, R1-R3 neighbour votes; clamps for hard flags, control libraries and the ubiquitous human ASV; dynamic denominator. | Design discussion. | not scored |
| v1 | 2026-10-03 | Clamp every detection of PCR positive-control component ASVs; cap U1's evidence for signal at +1 decade. | Gross positive-control contamination (e.g. *Lithobates* at high abundance in Brazil and Sierra Nevada) cannot be told from signal by read structure; absence from a few negatives is weak evidence. | not scored |
| v2 | before 2026-10-05 | Cap U4's evidence for signal (logit 1); add 2,000 pseudo reads to field-library denominators; make R3 one-sided (sink evidence only). | Absence from negatives is weak evidence; a high x resting on few remaining reads is weak evidence; the library vote was snowballing towards signal. | not scored |
| v2.1 | 2026-10-07 | Clamp the reviewed `known_contaminants` list instead of the ubiquitous human ASV. | Human DNA may be real signal in field water; lab stocks and positive-control constructs are contamination everywhere. | not scored |
| ladder-1 | 2026-10-08 | Model ladder B1-M5 fitted on the development labels as planned above (`prior_sd` 2.5); legacy moved to `13_legacy_edge_classification.Rmd`. | Planned before labelling. | not scored |

v0-v2.1 were made before this protocol: they were motivated partly by inspecting individual cases, some judged with GBIF locality. They are kept as the `legacy` baseline. The ladder above replaces them as the main line of development; its rungs are logged here as they are fitted.

## After Preliminary Analysis

You are now ready to move into your own analysis to address the research questions you have at hand!

**Note on archiving run directories:** Once you are confident that you are done with the processing and analysis workflows for a given dataset, you may consider deleting the `sequences` and `analysis` folders from your run directory. The remaining `outputs` folder should have all the outputs from theis pipeline, and metadata on the sequence files used and the configuration settings if needed for future reference. This will help you reclaim some space on your machine (just make sure that the sequence files are archived somewhere else!).
