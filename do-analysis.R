#!/usr/bin/env Rscript
## =============================================================================
## do-analysis.R -- run the full Wales 20 mph evaluation pipeline, in order.
##
## Usage:
##   Rscript do-analysis.R                 # run everything (see timings below)
##   Rscript do-analysis.R panel analysis  # run selected stages only
##   Rscript do-analysis.R list            # show available stages and inputs
##
## Stages (in dependency order):
##   prep      regenerate rectangular CSVs from the ODS downloads (python3)
##   extract   build collisions_2017_2025.csv / casualties_2017_2025.csv from
##             the raw STATS19 all-years downloads (requires awk via the shell)
##   panel     scripts/build_panel.R         -- data/panel.rds (~3 min)
##   analysis  scripts/analysis.R            -- headline DiD, event study, KSI
##                                             permutation (~75 min)
##   perm_all  scripts/permutation_all.R    -- all-casualty permutation (~1 min)
##   esfull    scripts/event_study_full.R   -- full-pre-period event studies for
##                                             Figures 1/1a/combined (~15 min)
##   itsann    scripts/its_annotations.R    -- ITS slopes/steps for fig1_its
##                                             annotations (~15 min)
##   modelfit  scripts/model_fit.R          -- fit metrics for fig3 labels (~15 min)
##   occupant  scripts/occupant_effects.R   -- vehicle-occupant DiD effects (~5 min)
##   mech      scripts/mechanism_decomposition.R -- collision/severity chain (~15 min)
##   dose      scripts/dose_response.R       -- treatment-intensity scaling (~10 min)
##   its       scripts/its_decomposition.R    -- level/slope decomposition (~20 min)
##   map       scripts/fix_map_model.R       -- LA-level adjusted-KSI effects (~8 min)
##   mc        scripts/multiple_comparisons.R-- Holm + BH corrections (seconds)
##   precovid  scripts/robustness_precovid.R -- 2017-19 vs 2024 spec (~15 min)
##   figures   scripts/figures.R             -- all figures (~3 min)
##
## Prerequisites (download locations in README.md "Sources"):
##   data/raw/stats19/collision_all_years.csv   (DfT STATS19, 1979-latest)
##   data/raw/stats19/casualty_all_years.csv     (DfT STATS19, 1979-latest)
##   data/raw/traffic/tra8904-km-by-la.ods       (DfT TRA8904)
##   data/raw/traffic/tra0413-km-pedal-cycle-traffic-by-region.ods (DfT TRA0413)
##   data/raw/pop/myebtablesuk20112025.xlsx     (ONS MYEB2)
##   data/raw/lad_to_county_apr2023.csv          (ONS Open Geography)
##   data/raw/geo/lad_boundaries.geojson        (ONS Open Geography)
##
## Each stage is run as a separate Rscript process; the wrapper stops at the
## first failure and reports which stage failed and how long each stage took.
## =============================================================================

stages <- c("prep", "extract", "panel", "analysis", "perm_all", "esfull",
            "itsann", "modelfit", "occupant", "its",
            "map", "mc", "mech", "dose", "precovid", "figures")

## ---- argument parsing -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- "all"
if (identical(args, "list")) {
  cat("Stages:\n")
  for (s in stages) cat(" ", s, "\n", sep = "")
  cat("  all (default) = all stages in order\n")
  quit(save = "no", status = 0)
}
if (identical(args, "all")) args <- stages
bad <- setdiff(args, stages)
if (length(bad) > 0) stop("Unknown stage(s): ", paste(bad, collapse = ", "),
                         "\nRun Rscript do-analysis.R list")
args <- stages[stages %in% args]   # preserve dependency order regardless

## ---- helpers ----------------------------------------------------------------
required_raw <- c(
  "data/raw/stats19/collision_all_years.csv",
  "data/raw/stats19/casualty_all_years.csv",
  "data/raw/traffic/tra8904-km-by-la.ods",
  "data/raw/traffic/tra0413-km-pedal-cycle-traffic-by-region.ods",
  "data/raw/pop/myebtablesuk20112025.xlsx",
  "data/raw/lad_to_county_apr2023.csv",
  "data/raw/geo/lad_boundaries.geojson")

run_r <- function(script) {
  cat("\n=== Running", script, "===\n")
  t0 <- Sys.time()
  status <- suppressWarnings(system2("Rscript", shQuote(script),
                                      stdout = "", stderr = ""))
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  if (!identical(status, 0)) {
    stop("Stage failed (exit status ", status, ") in ", script,
         " after ", dt, " min. Fix and re-run; completed stages are kept.")
  }
  cat("---", script, "done in", dt, "min ---\n")
}

run_shell <- function(cmd) {
  cat("Running:", cmd, "\n")
  status <- suppressWarnings(system(cmd))
  if (!identical(status, 0)) stop("Command failed (status ", status, "): ", cmd)
}

## ---- stage: prep (ODS -> rectangular CSV) -----------------------------------
do_prep <- function() {
  if (!file.exists("data/raw/traffic/tra8904_full.csv"))
    run_shell("python3 scripts/tra8904_to_csv.py")   # needs the TRA8904 ODS
  if (!file.exists("data/raw/traffic/tra0413_raw.csv"))
    run_shell(paste0("python3 scripts/ods_to_csv.py ",
                     "data/raw/traffic/tra0413-km-pedal-cycle-traffic-by-region.ods ",
                     "data/raw/traffic/tra0413_raw.csv"))
}

## ---- stage: extract (STATS19 all-years -> 2017-2025 analysis files) ---------
## These are the awk extractions used to build the panel inputs; they keep the
## panel build fast. If you prefer not to use the shell, point
## scripts/build_panel.R at your own extraction.
do_extract <- function() {
  if (file.exists("data/raw/stats19/collisions_2017_2025.csv") &&
      file.exists("data/raw/stats19/casualties_2017_2025.csv")) {
    cat("Extracted STATS19 files already present; skipping.\n")
    return(invisible())
  }
  run_shell(paste(
    "awk -F',' 'NR==1{print \"collision_index,collision_year,date,la_ons,road_type,first_road_class,speed_limit,urban_or_rural\"}",
    "$2>=2017 && $2<=2025 {print $1\",\"$2\",\"$12\",\"$16\",\"$21\",\"$19\",\"$22\",\"$37}'",
    "data/raw/stats19/collision_all_years.csv > data/raw/stats19/collisions_2017_2025.csv"))
  run_shell(paste(
    "awk -F',' 'NR==1{print \"collision_index,collision_year,casualty_class,age_of_casualty,casualty_severity,casualty_type,injury_based,adj_serious,adj_slight\"}",
    "$2>=2017 && $2<=2025 {print $1\",\"$2\",\"$6\",\"$8\",\"$10\",\"$16\",\"$20\",\"$21\",\"$22}'",
    "data/raw/stats19/casualty_all_years.csv > data/raw/stats19/casualties_2017_2025.csv"))
}

## ---- stage runners ----------------------------------------------------------
stage_fns <- list(
  prep     = do_prep,
  extract  = do_extract,
  panel    = function() run_r("scripts/build_panel.R"),
  analysis = function() run_r("scripts/analysis.R"),
  perm_all = function() run_r("scripts/permutation_all.R"),
  esfull   = function() run_r("scripts/event_study_full.R"),
  itsann   = function() run_r("scripts/its_annotations.R"),
  modelfit = function() run_r("scripts/model_fit.R"),
  occupant = function() run_r("scripts/occupant_effects.R"),
  mech     = function() run_r("scripts/mechanism_decomposition.R"),
  dose     = function() run_r("scripts/dose_response.R"),
  its      = function() run_r("scripts/its_decomposition.R"),
  map      = function() run_r("scripts/fix_map_model.R"),
  mc       = function() run_r("scripts/multiple_comparisons.R"),
  precovid = function() run_r("scripts/robustness_precovid.R"),
  figures  = function() run_r("scripts/figures.R"))

## ---- pre-flight: raw downloads must exist ------------------------------------
missing <- required_raw[!file.exists(required_raw)]
if (length(missing) > 0)
  stop("Missing raw inputs (see README.md Sources):\n  ",
       paste(missing, collapse = "\n  "))

## derived inputs: add prep/extract stages if their outputs are missing
need_prep <- !file.exists("data/raw/traffic/tra8904_full.csv") ||
            !file.exists("data/raw/traffic/tra0413_raw.csv")
if (need_prep && !("prep" %in% args)) {
  cat("Note: derived traffic CSVs missing; adding 'prep' stage.\n")
  args <- stages[stages %in% c(args, "prep")]
}
need_extract <- !file.exists("data/raw/stats19/collisions_2017_2025.csv") ||
                !file.exists("data/raw/stats19/casualties_2017_2025.csv")
if (need_extract && !("extract" %in% args)) {
  cat("Note: extracted STATS19 files missing; adding 'extract' stage.\n")
  args <- stages[stages %in% c(args, "extract")]
}

## ---- run --------------------------------------------------------------------
cat("Wales 20 mph evaluation pipeline\n")
cat("Stages requested:", paste(args, collapse = " -> "), "\n")
if (all(c("analysis", "perm_all", "esfull", "itsann", "modelfit", "occupant",
          "its", "map", "mc", "mech", "dose", "precovid") %in% args)) {
  cat("Note: full run takes roughly 2.5 hours of compute.\n")
}
t_total <- Sys.time()
for (s in args) stage_fns[[s]]()
cat("\nAll requested stages completed in",
    round(as.numeric(difftime(Sys.time(), t_total, units = "mins")), 1),
    "min.\n")
