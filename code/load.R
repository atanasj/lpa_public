### load.R: Takes care of loading in all the data required. Typically this is a
### short file, reading in data from files, URLs and/or ODBC. Depending on the
### project at this point I'll either write out the workspace using save() or
### just keep things in memory for the next step.

### ============================================================================
### LOAD DATA AND SET OUTPUT NAMES
### ============================================================================

## NOTE all locations are relative to this project for use with the `here` lib
## data objects to load here
## get original dataset
## NOTE this data set is corrected with added reverse score vars (suffix 'r')
.phd_sp_full_df <-
  Sys.glob("data/raw/sp_long_full_fix-dup_2021*.csv")

## data set of participants with at least first session
.phd_sp_df <-
  "data/proc/phd_sp_v2.rds"

## data set used to fit the lpa
.lpa_fit_df <-
  "data/proc/lpa_fit_v2.rds"

## data set used to fit the lta
.lta_fit_df <-
  "data/proc/lta_fit.rds"

## list of data for optseed specification stage
.m1_optseed <-
  "data/proc/m1_optseed.rds"

.m1_optseed_rv <-
  "data/proc/m1_optseed_rv.rds"

## list of data for LMR and BLRT specification stage
.m2_ll_tests <-
  "data/proc/m2_ll_tests.rds"

.m2_ll_nohx <-
  "data/proc/m2_ll_nohx.rds"

.m2_ll_hisi <-
  "data/proc/m2_ll_hisi.rds"

.m2_ll_si_nosa <-
  "data/proc/m2_ll_si_nosa.rds"

.m2_ll_si_sa <-
  "data/proc/m2_ll_si_sa.rds"

.m2_lpa_si <-
  "data/proc/m2_lpa_si.rds"

.m2_lpa_si2 <-
  "data/proc/m2_lpa_si2.rds"

.m2_ll_tests_rv <-
  "data/proc/m2_ll_tests_rv.rds"

.lpa_fin <-
  "data/proc/lpa_fin.rds"

.m2_lta_t0 <-
  "data/proc/m2_lta_t0.rds"

.m2_lta_t1 <-
  "data/proc/m2_lta_t1.rds"

### ==========================================================================
### END SCRIPT: load.R
### ============================================================================
