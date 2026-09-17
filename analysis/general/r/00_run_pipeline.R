# Runs the classification/export pipeline (scripts 04-09) end to end for one run.
# Edit pipeline_config.yml, then source this script. Each script reads pipeline_config.yml
# itself, so the same config applies whether a script is run manually (open it in RStudio and
# run chunks) or via this driver -- there are no separate params to keep in sync.
# Script 05 (GBIF query) only runs when pipeline_config.yml sets gbif_query: true.
#
# Each script's own write_csv()/write.csv() calls save the actual pipeline outputs to
# run_dir/output as usual. Knitting still produces a throwaway HTML doc per script (that's
# how rmarkdown::render executes the chunks) -- those are written to a session temp dir and
# are not meant to be kept.

librarian::shelf(here, yaml, rmarkdown)

script_dir = here("analysis", "general", "r")
run_config = read_yaml(here(script_dir, "pipeline_config.yml"))

scratch_dir = file.path(tempdir(), "edna_pipeline_run")
dir.create(scratch_dir, showWarnings = FALSE, recursive = TRUE)

render_step = function(rmd) {
  cat("\n==== Running", rmd, "====\n")
  rmarkdown::render(
    input = file.path(script_dir, rmd),
    envir = new.env(),
    knit_root_dir = here(),
    output_dir = scratch_dir,
    quiet = TRUE
  )
}

render_step("04_web_blast_json_parse.Rmd")

if (isTRUE(run_config$gbif_query)) {
  render_step("05_query_taxonomy_geography.Rmd")
}

render_step("06_classify_asv.Rmd")
render_step("07_sample_map.Rmd")
render_step("08_sample_controls.Rmd")
render_step("09_export_results.Rmd")

cat("\nPipeline complete. Outputs saved under:", here(run_config$run_dir, "output"), "\n")
