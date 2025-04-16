### ============================================================================
### SOURCE SCRIPTS
### ============================================================================
rm(list = ls(all = TRUE))
source("./code/func.R")
source("./code/load.R")

### ============================================================================
### LOAD DATA
### ============================================================================

.lta_fit <-
  readRDS(here(.lta_fit_df))

lta_fit_0 <-
  .lta_fit %>%
  filter(session == 1) %>%
  select(ID, pb, tb, bhs, ac)

lta_fit_1 <-
  .lta_fit %>%
  filter(final_session == 1) %>%
  select(ID, pb, tb, bhs, ac)

## load data with revers scroed vars as well as original vars
.lpa_fit <-
  readRDS(here(.lpa_fit_df)) %>%
  mutate(
    hxsa = case_when(
      suicide_attempts <= 0 ~ 0,
      suicide_attempts >= 1 ~ 1,
      TRUE ~ NA_real_
    )
  )

lpa_fit <-
  .lpa_fit %>%
  select(ID, pb, tb, bhs, ac)

lpa_fit_nohxsa <-
  .lpa_fit %>%
  filter(hxsa %in% 0) %>%
  select(ID, pb, tb, bhs, ac)

lpa_fit_hisi <-
  .lpa_fit %>%
  filter(si > 8) %>%
  select(ID, pb, tb, bhs, ac)

lpa_fit_si_nosa <-
  .lpa_fit %>%
  filter(si > 8 & hxsa %in% 0) %>%
  select(ID, pb, tb, bhs, ac)

lpa_fit_si_sa <-
  .lpa_fit %>%
  filter(si > 8 & hxsa %in% 1) %>%
  select(ID, pb, tb, bhs, ac)

lpa_si <-
  .lpa_fit %>%
  filter(session == 1, final_session != 1) %>%
  filter(si > 0 | siss > 0) %>%
  select(ID, pb, tb, bhs, ac)

lpa_si2 <-
  .lpa_fit %>%
  filter(session == 1, final_session != 1) %>%
  filter(si > 0) %>%
  select(ID, pb, tb, bhs, ac)

## load full data for to merge with profiles
phd_sp_fit <-
  readRDS(here(.phd_sp_df))

phd_sp_fit %>%
  group_by(ID) %>%
  mutate(
    max_ses = max(session, na.rm = TRUE),
    pf = max(planned_final, na.rm = TRUE),
  ) %>%
  ungroup() %>%
  mutate(
    nsi = case_when(
      si == 0 ~ 0,
      si > 0  ~ 1,
      TRUE ~ NA_real_
    ),
    nsi2 = case_when(
      siss == 0 ~ 0,
      siss > 0  ~ 1,
      TRUE ~ NA_real_
    ),
    nsi3 = case_when(
      hxsuicide == 0 ~ 0,
      hxsuicide > 0  ~ 1,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(session == 1, final_session != 1) %>%
  rowwise() %>%
  mutate(
    nsi4 = if_else(sum(c_across(nsi:nsi2), na.rm = TRUE) > 0, 1, 0)
  ) %>%
  ## group_by(nsi, nsi2, nsi3) %>%
  group_by(nsi4) %>%
  summarise(
    n = n(),
    across(
      c(si, siss, max_ses, pf),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(.x, na.rm = TRUE)
      ),
      .names = "{.col}.{.fn}"
    )
  ) %>%
  pivot_longer(
    ## !nsi2
    !nsi4
  ) %>%
  pivot_wider(
    ## names_from = nsi2
    names_from = nsi4
  ) %>%
  print(n = Inf)



### ============================================================================
### NOTES ON M1_OPTSEED
### ============================================================================

## NOTE the process for specifying the models is:
## =============================================================================
## 01. use function below to create an initial `.inp` file in its own directory
##     without running the TECH11 or TECH14 commands
## 02. run each `.inp` file via the cli e.g., `mplus c{k,modelname}.inp
##     &>/dev/null`
## 03. modify each `.inp` file increasing the `STARTS` value until the best
##     loglikelihood s been repeated
## 04. repeat this with each `.inp` until statisfied above
## 05. copy with above
## 06. set `OPTSEED` value for each `inp` file
## 07. run initially via mplusAutomation to get the initial `.out` files
## 08. update the `LRTSTARTS` value as per the `.out` file recommendations per
##     `.inp` file
## 09. run the readModels function on the directory
## 10. once the number of latent classes are estimated, then merge the outputted
##     data of the best fitting model and merge back with the original data set.
## 11. run each follow up analysis within R.

## LPA STEP 0.1: MISSING DATA: FULL SAMPLE
## =============================================================================
lpa_miss <-
  mplusObject(
    TITLE = "LPA missing data",
    VARIABLE = glue(
      "USEVARIABLES = pb tb bhs ac;
       IDVARIABLE = ID;"
    ),
    ANALYSIS =
      "ESTIMATOR = MLR;
       ALGORITHM = INTEGRATION;
       TYPE = BASIC;
       ! STARTS = 160 32;
       ! OPTSEED = 391179;
       ! Dual-Core with hyperthreading = 4 processors
       PROCESSORS = 4;
       ! K-1STARTS = 20 4; ! used with TECH11 output
       ! LRTSTARTS = 20 4 40 8; ! used with TECH14 output",
    OUTPUT =
      "TECH1 TECH8;
       ! TECH11 TECH14;
       STANDARDIZED; ! this saves standardised scores for plotting",
    SAVEDATA =
      c(
        glue(
          "FILE = ",
          here(
            "mplus/1_optseed_m",
            "c0_optseed_m_v1.dat;"
          )
        ),
        "SAVE = CPROB;
         FORMAT = FREE;"
      ),
    PLOT =
      "TYPE = PLOT3;",
    usevariables = colnames(lpa_fit),
    rdata = lpa_fit
  )

lpa_miss_fit <-
  mplusModeler(lpa_miss,
    dataout = glue(here(
      "mplus/1_optseed_m", "full_1_optseed_m.dat"
    )),
    modelout = glue(here(
      "mplus/1_optseed_m", "c0_optseed_m_v1.inp"
    )),
    check = TRUE, run = TRUE, hashfilename = FALSE
  )

## LPA STEP 1.1: DETERMINE OPTSEED_V1: FULL SAMPLE
## =============================================================================
## k2 to k6 models explored
lpa_k26_1 <-
  ## how many classes do we want to estimate?
  lapply(c(7), function(k) {
    lpa_enum_optseed_v1 <-
      mplusObject(
        TITLE = glue("Class {k} 1_optseed_m_v2"),
        VARIABLE = glue(
          "USEVARIABLES = pb tb bhs ac;
           IDVARIABLE = ID;
           CLASSES = c({k});"
        ),
        ANALYSIS =
          "ESTIMATOR = MLR;
           ALGORITHM = INTEGRATION;
           TYPE = MIXTURE;
           STARTS = 80 16;
           ! OPTSEED = 391179;
           ! Dual-Core with hyperthreading = 4 processors
           PROCESSORS = 4;
           ! K-1STARTS = 20 4; ! used with TECH11 output
           ! LRTSTARTS = 20 4 40 8; ! used with TECH14 output",
        OUTPUT =
        ## TECH 1 = parameter spec & starting values for est parameter
        ## TECH 8 = provides the optimization history
        ## TECH 11 = Lo-Mendell-Rubin Adjusted LRT p-value
        ## TECH 14 = Bootstrapped Likelihood Ratio Test
        "TECH1 TECH8;
           ! TECH11 TECH14;
          STANDARDIZED; ! this saves standardised scores for plotting",
        SAVEDATA =
          c(
            glue(
              "FILE = ",
              here(
                "mplus/1_optseed_m",
                "c{k}_1_optseed_m_v2.dat;"
              )
            ),
            "SAVE = CPROB;
            FORMAT = FREE;"
          ),
        PLOT =
          "TYPE = PLOT3;",
        usevariables = colnames(lpa_fit),
        rdata = lpa_fit
      )
    lpa_miss_fit <-
      mplusModeler(lpa_enum_optseed_v1,
        dataout = glue(here(
          "mplus/1_optseed_m", "full_1_optseed_m.dat"
        )),
        modelout = glue(here(
          "mplus/1_optseed_m", "c{k}_1_optseed_m_v2.inp"
        )),
        check = TRUE, run = TRUE, hashfilename = FALSE
      )
  })

## LPA STEP 1.2: DETERMINE OPTSEED_V1: NO HXSA
## =============================================================================
## k2 to k6 models explored
lpa_k26_1 <-
  ## how many classes do we want to estimate?
  lapply(c(1:4, 6:7), function(k) {
    lpa_enum_optseed_v1 <-
      mplusObject(
        TITLE = glue("Class {k} 1_optseed_nohxsa_v2"),
        VARIABLE = glue(
          "USEVARIABLES = pb tb bhs ac;
           IDVARIABLE = ID;
           CLASSES = c({k});"
        ),
        ANALYSIS =
          "ESTIMATOR = MLR;
           ALGORITHM = INTEGRATION;
           TYPE = MIXTURE;
           STARTS = 40 8;
           ! OPTSEED = 391179;
           ! Dual-Core with hyperthreading = 4 processors
           PROCESSORS = 4;
           ! K-1STARTS = 20 4; ! used with TECH11 output
           ! LRTSTARTS = 20 4 40 8; ! used with TECH14 output",
        OUTPUT =
        ## TECH 1 = parameter spec & starting values for est parameter
        ## TECH 8 = provides the optimization history, remove after optseed set
        ## TECH 11 = Lo-Mendell-Rubin Adjusted LRT p-value
        ## TECH 14 = Bootstrapped Likelihood Ratio Test
        "TECH1 TECH8;
         ! TECH11 TECH14;
         STANDARDIZED; ! this saves standardised scores for plotting",
        SAVEDATA =
          c(
            glue(
              "FILE = ",
              here(
                "mplus/1_optseed_nohxsa",
                "c{k}_1_optseed_nohxsa_v2.dat;"
              )
            ),
            "SAVE = CPROB;
            FORMAT = FREE;"
          ),
        PLOT =
          "TYPE = PLOT3;",
        usevariables = colnames(lpa_fit_nohxsa),
        rdata = lpa_fit_nohxsa
      )
    lpa_miss_fit <-
      mplusModeler(lpa_enum_optseed_v1,
        dataout = glue(here(
          "mplus/1_optseed_nohxsa", "full_1_optseed_nohxsa.dat"
        )),
        modelout = glue(here(
          "mplus/1_optseed_nohxsa", "c{k}_1_optseed_nohxsa_v2.inp"
        )),
        check = TRUE, run = TRUE, hashfilename = FALSE
      )
  })

## run all models with optseed added
runModels("./mplus/2_ll_nohx")

## LPA STEP 1.3: DETERMINE OPTSEED_V1: HISI
## =============================================================================
## k2 to k6 models explored
lpa_k26_1 <-
  ## how many classes do we want to estimate?
  lapply(c(5), function(k) {
    lpa_enum_optseed_v1 <-
      mplusObject(
        TITLE = glue("Class {k} 1_optseed_hisi_v2"),
        VARIABLE = glue(
          "USEVARIABLES = pb tb bhs ac;
           IDVARIABLE = ID;
           CLASSES = c({k});"
        ),
        ANALYSIS =
          "ESTIMATOR = MLR;
           ALGORITHM = INTEGRATION;
           TYPE = MIXTURE;
           STARTS = 160 32;
           ! OPTSEED = 391179;
           ! Dual-Core with hyperthreading = 4 processors
           PROCESSORS = 4;
           ! K-1STARTS = 20 4; ! used with TECH11 output
           ! LRTSTARTS = 20 4 40 8; ! used with TECH14 output",
        OUTPUT =
        ## TECH 1 = parameter spec & starting values for est parameter
        ## TECH 8 = provides the optimization history, remove after optseed set
        ## TECH 11 = Lo-Mendell-Rubin Adjusted LRT p-value
        ## TECH 14 = Bootstrapped Likelihood Ratio Test
        "TECH1 TECH8;
         ! TECH11 TECH14;
         STANDARDIZED; ! this saves standardised scores for plotting",
        SAVEDATA =
          c(
            glue(
              "FILE = ",
              here(
                "mplus/1_optseed_hisi",
                "c{k}_1_optseed_hisi_v2.dat;"
              )
            ),
            "SAVE = CPROB;
            FORMAT = FREE;"
          ),
        PLOT =
          "TYPE = PLOT3;",
        usevariables = colnames(lpa_fit_hisi),
        rdata = lpa_fit_hisi
      )
    lpa_miss_fit <-
      mplusModeler(lpa_enum_optseed_v1,
        dataout = glue(here(
          "mplus/1_optseed_hisi", "full_1_optseed_hisi.dat"
        )),
        modelout = glue(here(
          "mplus/1_optseed_hisi", "c{k}_1_optseed_hisi_v2.inp"
        )),
        check = TRUE, run = TRUE, hashfilename = FALSE
      )
  })

## run all models with optseed added
runModels("./mplus/2_ll_hisi")

## LPA STEP 1.4: DETERMINE OPTSEED_V1: HISI NO SA
## =============================================================================
## k2 to k6 models explored
lpa_k26_1 <-
  ## how many classes do we want to estimate?
  lapply(c(7), function(k) {
    lpa_enum_optseed_v1 <-
      mplusObject(
        TITLE = glue("Class {k} 1_opt_si_nosa_v2"),
        VARIABLE = glue(
          "USEVARIABLES = pb tb bhs ac;
           IDVARIABLE = ID;
           CLASSES = c({k});"
        ),
        ANALYSIS =
          "ESTIMATOR = MLR;
           ALGORITHM = INTEGRATION;
           TYPE = MIXTURE;
           STARTS = 160 64;
           ! OPTSEED = 391179;
           ! Dual-Core with hyperthreading = 4 processors
           PROCESSORS = 4;
           ! K-1STARTS = 20 4; ! used with TECH11 output
           ! LRTSTARTS = 20 4 40 8; ! used with TECH14 output",
        OUTPUT =
        ## TECH 1 = parameter spec & starting values for est parameter
        ## TECH 8 = provides the optimization history, remove after optseed set
        ## TECH 11 = Lo-Mendell-Rubin Adjusted LRT p-value
        ## TECH 14 = Bootstrapped Likelihood Ratio Test
        "TECH1 TECH8;
         ! TECH11 TECH14;
         STANDARDIZED; ! this saves standardised scores for plotting",
        SAVEDATA =
          c(
            glue(
              "FILE = ",
              here(
                "mplus/1_opt_si_nosa",
                "c{k}_1_opt_si_nosa_v2.dat;"
              )
            ),
            "SAVE = CPROB;
            FORMAT = FREE;"
          ),
        PLOT =
          "TYPE = PLOT3;",
        usevariables = colnames(lpa_fit_si_nosa),
        rdata = lpa_fit_si_nosa
      )
    lpa_miss_fit <-
      mplusModeler(lpa_enum_optseed_v1,
        dataout = glue(here(
          "mplus/1_opt_si_nosa", "full_1_opt_si_nosa.dat"
        )),
        modelout = glue(here(
          "mplus/1_opt_si_nosa", "c{k}_1_opt_si_nosa_v2.inp"
        )),
        check = TRUE, run = TRUE, hashfilename = FALSE
      )
  })

## run all models with optseed added
runModels("./mplus/2_ll_si_nosa")

## LPA STEP 1.5: DETERMINE OPTSEED_V1: HISI SA
## =============================================================================
## k2 to k6 models explored
lpa_k26_1 <-
  ## how many classes do we want to estimate?
  lapply(c(5,7), function(k) {
    lpa_enum_optseed_v1 <-
      mplusObject(
        TITLE = glue("Class {k} 1_opt_si_sa_v2"),
        VARIABLE = glue(
          "USEVARIABLES = pb tb bhs ac;
           IDVARIABLE = ID;
           CLASSES = c({k});"
        ),
        ANALYSIS =
          "ESTIMATOR = MLR;
           ALGORITHM = INTEGRATION;
           TYPE = MIXTURE;
           STARTS = 160 32;
           ! OPTSEED = 391179;
           ! Dual-Core with hyperthreading = 4 processors
           PROCESSORS = 4;
           ! K-1STARTS = 20 4; ! used with TECH11 output
           ! LRTSTARTS = 20 4 40 8; ! used with TECH14 output",
        OUTPUT =
        ## TECH 1 = parameter spec & starting values for est parameter
        ## TECH 8 = provides the optimization history, remove after optseed set
        ## TECH 11 = Lo-Mendell-Rubin Adjusted LRT p-value
        ## TECH 14 = Bootstrapped Likelihood Ratio Test
        "TECH1 TECH8;
         ! TECH11 TECH14;
         STANDARDIZED; ! this saves standardised scores for plotting",
        SAVEDATA =
          c(
            glue(
              "FILE = ",
              here(
                "mplus/1_opt_si_sa",
                "c{k}_1_opt_si_sa_v2.dat;"
              )
            ),
            "SAVE = CPROB;
            FORMAT = FREE;"
          ),
        PLOT =
          "TYPE = PLOT3;",
        usevariables = colnames(lpa_fit_si_sa),
        rdata = lpa_fit_si_sa
      )
    lpa_miss_fit <-
      mplusModeler(lpa_enum_optseed_v1,
        dataout = glue(here(
          "mplus/1_opt_si_sa", "full_1_opt_si_sa.dat"
        )),
        modelout = glue(here(
          "mplus/1_opt_si_sa", "c{k}_1_opt_si_sa_v2.inp"
        )),
        check = TRUE, run = TRUE, hashfilename = FALSE
      )
  })

## run all models with optseed added
runModels("./mplus/2_ll_si_sa")

## LPA STEP 1.1: DETERMINE OPTSEED_V1: SI NON-ZERO
## =============================================================================
## k2 to k6 models explored
lpa_k26_1 <-
  ## how many classes do we want to estimate?
  lapply(c(2:6), function(k) {
    lpa_enum_optseed_v1 <-
      mplusObject(
        TITLE = glue("Class {k} lpa_si_v2"),
        VARIABLE = glue(
          "USEVARIABLES = pb tb bhs ac;
           IDVARIABLE = ID;
           CLASSES = c({k});"
        ),
        ANALYSIS =
          "ESTIMATOR = MLR;
           ALGORITHM = INTEGRATION;
           TYPE = MIXTURE;
           STARTS = 40 8;
           ! OPTSEED = 391179;
           ! Dual-Core with hyperthreading = 4 processors
           PROCESSORS = 4;
           ! K-1STARTS = 20 4; ! used with TECH11 output
           ! LRTSTARTS = 20 4 40 8; ! used with TECH14 output",
        OUTPUT =
        ## TECH 1 = parameter spec & starting values for est parameter
        ## TECH 8 = provides the optimization history
        ## TECH 11 = Lo-Mendell-Rubin Adjusted LRT p-value
        ## TECH 14 = Bootstrapped Likelihood Ratio Test
          "TECH1; ! TECH8;
          ! TECH11 TECH14;
          STANDARDIZED; ! this saves standardised scores for plotting",
        SAVEDATA =
          c(
            glue(
              "FILE = ",
              here(
                "mplus/lpa_si",
                "c{k}_lpa_si_v2.dat;"
              )
            ),
            "SAVE = CPROB;
            FORMAT = FREE;"
          ),
        PLOT =
          "TYPE = PLOT3;",
        usevariables = colnames(lpa_si),
        rdata = lpa_si
      )
    lpa_miss_fit <-
      mplusModeler(lpa_enum_optseed_v1,
        dataout = glue(here(
          "mplus/lpa_si", "full_lpa_si.dat"
        )),
        modelout = glue(here(
          "mplus/lpa_si", "c{k}_lpa_si_v2.inp"
        )),
        check = TRUE, run = TRUE, hashfilename = FALSE
      )
  })

runModels("./mplus/2_lpa_si")


## LPA STEP 1.1: DETERMINE OPTSEED_V1: MSSI NON-ZERO
## =============================================================================
## k2 to k6 models explored
lpa_k26_1 <-
  ## how many classes do we want to estimate?
  lapply(c(1:4), function(k) {
    lpa_enum_optseed_v1 <-
      mplusObject(
        TITLE = glue("Class {k} lpa_si2_v2"),
        VARIABLE = glue(
          "USEVARIABLES = pb tb bhs ac;
           IDVARIABLE = ID;
           CLASSES = c({k});"
        ),
        ANALYSIS =
          "ESTIMATOR = MLR;
           ALGORITHM = INTEGRATION;
           TYPE = MIXTURE;
           STARTS = 40 8;
           ! OPTSEED = 391179;
           ! Dual-Core with hyperthreading = 4 processors
           PROCESSORS = 4;
           ! K-1STARTS = 20 4; ! used with TECH11 output
           ! LRTSTARTS = 20 4 40 8; ! used with TECH14 output",
        OUTPUT =
        ## TECH 1 = parameter spec & starting values for est parameter
        ## TECH 8 = provides the optimization history
        ## TECH 11 = Lo-Mendell-Rubin Adjusted LRT p-value
        ## TECH 14 = Bootstrapped Likelihood Ratio Test
          "TECH1; ! TECH8;
          ! TECH11 TECH14;
          STANDARDIZED; ! this saves standardised scores for plotting",
        SAVEDATA =
          c(
            glue(
              "FILE = ",
              here(
                "mplus/lpa_si2",
                "c{k}_lpa_si2_v2.dat;"
              )
            ),
            "SAVE = CPROB;
            FORMAT = FREE;"
          ),
        PLOT =
          "TYPE = PLOT3;",
        usevariables = colnames(lpa_si2),
        rdata = lpa_si2
      )
    lpa_miss_fit <-
      mplusModeler(lpa_enum_optseed_v1,
        dataout = glue(here(
          "mplus/lpa_si2", "full_lpa_si2.dat"
        )),
        modelout = glue(here(
          "mplus/lpa_si2", "c{k}_lpa_si2_v2.inp"
        )),
        check = TRUE, run = TRUE, hashfilename = FALSE
      )
  })

runModels("./mplus/2_lpa_si2")

## READ IN ALL THE MODELS (V1 AND V2)
## -----------------------------------------------------------------------------
optseed_m <-
  readModels(here("mplus/1_optseed_m"), quiet = TRUE)

optseed_m %>%
  saveRDS(here(.m1_optseed))

## READ IN ALL THE MODELS WITH TECH11 TECH14 OUTPUT
## -----------------------------------------------------------------------------
ll_tests_m <-
  readModels(here("mplus/2_ll_tests_m"), quiet = TRUE)

ll_tests_m %>%
  saveRDS(here(.m2_ll_tests))

ll_nohx <-
  readModels(here("mplus/2_ll_nohx"), quiet = TRUE)

ll_nohx %>%
  saveRDS(here(.m2_ll_nohx))

ll_hisi <-
  readModels(here("mplus/2_ll_hisi"), quiet = TRUE)

ll_hisi %>%
  saveRDS(here(.m2_ll_hisi))

ll_si_nosa <-
  readModels(here("mplus/2_ll_si_nosa"), quiet = TRUE)

ll_si_nosa %>%
  saveRDS(here(.m2_ll_si_nosa))

ll_si_sa <-
  readModels(here("mplus/2_ll_si_sa"), quiet = TRUE)

ll_si_sa %>%
  saveRDS(here(.m2_ll_si_sa))

lpa_si <-
  readModels(here("mplus/2_lpa_si"), quiet = TRUE)

lpa_si %>%
  saveRDS(here(.m2_lpa_si))


lpa_si2 <-
  readModels(here("mplus/2_lpa_si2"), quiet = TRUE)

lpa_si2 %>%
  saveRDS(here(.m2_lpa_si2))

### ============================================================================
### MODEL SELECTION
### ============================================================================
.lpa_class_select <-
  list(
    full = ll_tests_m,
    no_sa = ll_nohx,
    hi_si = ll_hisi,
    hi_si_no_sa = ll_si_nosa,
    hi_si_sa = ll_si_sa,
    si = lpa_si,
    si2 = lpa_si2
  ) %>%
  map(
    ~ .x %>%
      SummaryTable(
        keepCols = c(
          "Title",
          "LL",
          "Entropy",
          "BIC",        # most trusted stat
          "aBIC",       # remove from kbl
          "T11_KM1LL",
          "T11_VLMR_PValue",
          "BLRT_KM1LL", # remove from kbl
          "BLRT_PValue" # remove from kbl
        ),
        sortBy = "Title"
      ) %>%
      mutate(
        Classes = str_c(str_extract(Title, "\\d"), " Classes"), .before = 1,
        across(where(is.numeric) & !ends_with("Value"), ~ round(.x, 2)),
        across(ends_with("Value"), ~ round(.x, 3)),
        ## number of classes
        ## K = 2:4,
        ## BAYES FACTOR comparision of models with k vs k+1 classes
        BF = signif(exp(-0.5 * BIC - -0.5 * lead(BIC)), 2),
        ## SCHWARTZ INFOMRATION CRITERION (SIC)
        SIC = -0.5 * BIC,
        cmP = exp(SIC - max(SIC))
      ) %>%
      select(-Title, -SIC) %>%
      mutate(
        BF_fx = case_when(
          BF < 3 ~ "weak",
          BF > 3 & BF < 10 ~ "mod",
          BF > 10 ~ "strong",
          TRUE ~ NA_character_
        ),
        cmP_fx = if_else(cmP == max(cmP), "pref", "no_pref")
      ) %>%
      relocate(BF, .before = BF_fx) %>%
      relocate(cmP, .before = cmP_fx) %>%
      filter(!Classes %in% c("1 Classes", "7 Classes"))
  ) %>%
  bind_rows(., .id = "Sample")

.lpa_class_select_N <-
  list(
    full = ll_tests_m,
    no_sa = ll_nohx,
    hi_si = ll_hisi,
    hi_si_no_sa = ll_si_nosa,
    hi_si_sa = ll_si_sa,
    si = lpa_si,
    si2 = lpa_si2
  ) %>%
  map(
    `[`, (2)
  ) %>%
  unlist(recursive = FALSE) %>%
  map(`[`, "savedata") %>%
  map(
    ~ .x %>%
      as.data.frame() %>%
      nrow()
  ) %>%
  bind_rows() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(
    Sample = 1,
    N = 2
  ) %>%
  separate(Sample, c("Sample", "NOPA"), sep = "\\.") %>%
  select(-NOPA)

lpa_class_select <-
  .lpa_class_select %>%
  full_join(.lpa_class_select_N) %>%
  relocate(N, .after = Sample)

lpa_class_plots <-
  list(
    full = ll_tests_m,
    no_sa = ll_nohx,
    hi_si = ll_hisi,
    hi_si_no_sa = ll_si_nosa,
    hi_si_sa = ll_si_sa,
    si = lpa_si,
    si2 = lpa_si2
  ) %>%
  map(`[`, (2:5)) %>%
  map(
    ~ .x %>%
      plotMixtures(
        coefficients = "stdyx.standardized")
  ) %>%
  wrap_plots()

ggsave(
  "./figures/lpa_class_plots.png",
  lpa_class_plots
)

### ============================================================================
### CREATE FINAL LPA DATASET
### ============================================================================

.lpa_full <-
  ll_tests_m$c4_2_ll_tests_m_v1.out$savedata %>%
  clean_names() %>%
  rename(
    ID = id,
    c_full = c,
    c_fullprob1 = cprob1,
    c_fullprob2 = cprob2,
    c_fullprob3 = cprob3,
    c_fullprob4 = cprob4
    ) %>%
  select(ID, c_full, c_fullprob1:c_fullprob4)

.lpa_hisi <-
  ll_hisi$c2_2_ll_hisi_v1.out$savedata %>%
  clean_names() %>%
  rename(
    ID = id,
    c_hisi =  c,
    c_hisiprob1 = cprob1,
    c_hisiprob2 = cprob2
    ) %>%
  select(ID, c_hisi, c_hisiprob1:c_hisiprob2)

.lpa_si <-
  lpa_si2$c3_2_lpa_si2_v1.out$savedata %>%
  clean_names() %>%
  rename(
    ID = id
    ) %>%
  select(ID, c, cprob1:cprob3)

lpa_fin <-
  phd_sp_fit %>%
  full_join(.lpa_si) %>%
  full_join(.lpa_full) %>%
  full_join(.lpa_hisi) %>%
  select(
    ID:planned_final, pb:ac, dep:isas_social,
    atsi:suburb, everything()
  ) %>%
  relocate(starts_with("c_"), .after = ID) %>%
  relocate(c(c, cprob1:cprob3), .after = ID)

lpa_fin %>%
  saveRDS(here(.lpa_fin))

lpa_fin %>%
  filter(
    if_any(c(session, final_session), ~ . %in% 1),
    if_all(c(sp_episode, consent), ~ . %in% 1)
    ) %>%
  select(uci:final_session, contains("mssi")) %>%
  split(.$final_session) %>%
  map(
    ~ .x %>%
      group_by(mssi_scrn_miss2, mssi_co) %>%
      summarise(
        n = n()
      )
  )

### ============================================================================
### LPA: ALTERNATIVE IDENTIFICATION
### ============================================================================

lpa_fit_2 <-
  .lpa_fit %>%
  mutate(sex2 = if_else(sex %in% 2, NA_real_, sex), .after = sex) %>%
  select(
    ID,
    pb, tb, bhs, ac, si,
    dep, siss, anx, str,
    starts_with(c("hx", "past")),
    isas_3a, age, sex2
  )

## LTA STEP 1.1: DETERMINE OPTSEED_V1: TIME 0
## =============================================================================
## k2 to k6 models explored
lpa_k26_1 <-
  ## how many classes do we want to estimate?
  lapply(c(1:7), function(k) {
    lpa_enum_optseed_v1 <-
      mplusObject(
        TITLE = glue("Class {k} 1_lta_t0_v2"),
        VARIABLE = glue(
          "USEVARIABLES = pb, tb, bhs;
           IDVARIABLE = ID;
           CLASSES = c({k});"
        ),
        ANALYSIS =
          "ESTIMATOR = MLR;
           ALGORITHM = INTEGRATION;
           TYPE = MIXTURE;
           STARTS = 40 8;
           ! OPTSEED = 391179;
           ! Dual-Core with hyperthreading = 4 processors
           PROCESSORS = 4;
           ! K-1STARTS = 20 4; ! used with TECH11 output
           ! LRTSTARTS = 20 4 40 8; ! used with TECH14 output",
        OUTPUT =
        ## TECH 1 = parameter spec & starting values for est parameter
        ## TECH 8 = provides the optimization history, remove after optseed set
        ## TECH 11 = Lo-Mendell-Rubin Adjusted LRT p-value
        ## TECH 14 = Bootstrapped Likelihood Ratio Test
        "TECH1 TECH8;
         ! TECH11 TECH14;
         STANDARDIZED; ! this saves standardised scores for plotting",
        SAVEDATA =
          c(
            glue(
              "FILE = ",
              here(
                "mplus/1_lta_t0",
                "c{k}_1_lta_t0_v2.dat;"
              )
            ),
            "SAVE = CPROB;
            FORMAT = FREE;"
          ),
        PLOT =
          "TYPE = PLOT3;",
        usevariables = colnames(lta_fit_0),
        rdata = lta_fit_0
      )
    lpa_miss_fit <-
      mplusModeler(lpa_enum_optseed_v1,
        dataout = glue(here(
          "mplus/1_lta_t0", "full_1_lta_t0.dat"
        )),
        modelout = glue(here(
          "mplus/1_lta_t0", "c{k}_1_lta_t0_v2.inp"
        )),
        check = TRUE, run = TRUE, hashfilename = FALSE
      )
  })

runModels("./mplus/2_lta_t0")


## LTA STEP 1.2: DETERMINE OPTSEED_V1: TIME 1
## =============================================================================
## k2 to k6 models explored
lpa_k26_1 <-
  ## how many classes do we want to estimate?
  lapply(c(1:5, 7), function(k) {
    lpa_enum_optseed_v1 <-
      mplusObject(
        TITLE = glue("Class {k} 1_lta_t1_v2"),
        VARIABLE = glue(
          "USEVARIABLES = pb, tb, bhs;
           IDVARIABLE = ID;
           CLASSES = c({k});"
        ),
        ANALYSIS =
          "ESTIMATOR = MLR;
           ALGORITHM = INTEGRATION;
           TYPE = MIXTURE;
           STARTS = 40 8;
           ! OPTSEED = 391179;
           ! Dual-Core with hyperthreading = 4 processors
           PROCESSORS = 4;
           ! K-1STARTS = 20 4; ! used with TECH11 output
           ! LRTSTARTS = 20 4 40 8; ! used with TECH14 output",
        OUTPUT =
        ## TECH 1 = parameter spec & starting values for est parameter
        ## TECH 8 = provides the optimization history, remove after optseed set
        ## TECH 11 = Lo-Mendell-Rubin Adjusted LRT p-value
        ## TECH 14 = Bootstrapped Likelihood Ratio Test
        "TECH1 TECH8;
         ! TECH11 TECH14;
         STANDARDIZED; ! this saves standardised scores for plotting",
        SAVEDATA =
          c(
            glue(
              "FILE = ",
              here(
                "mplus/1_lta_t1",
                "c{k}_1_lta_t1_v2.dat;"
              )
            ),
            "SAVE = CPROB;
            FORMAT = FREE;"
          ),
        PLOT =
          "TYPE = PLOT3;",
        usevariables = colnames(lta_fit_1),
        rdata = lta_fit_1
      )
    lpa_miss_fit <-
      mplusModeler(lpa_enum_optseed_v1,
        dataout = glue(here(
          "mplus/1_lta_t1", "full_1_lta_t1.dat"
        )),
        modelout = glue(here(
          "mplus/1_lta_t1", "c{k}_1_lta_t1_v2.inp"
        )),
        check = TRUE, run = TRUE, hashfilename = FALSE
      )
  })

## run all models with optseed added
runModels("./mplus/2_lta_t1")

### LTA MODEL SAVING
### ============================================================================

lta_t0 <-
  readModels(here("mplus/2_lta_t0"), quiet = TRUE)

lta_t0 %>%
  saveRDS(here(.m2_lta_t0))

lta_t1 <-
  readModels(here("mplus/2_lta_t1"), quiet = TRUE)

lta_t1 %>%
  saveRDS(here(.m2_lta_t1))

### ============================================================================
### LTA MODEL SELECTION
### ============================================================================
.lta_class_select <-
  list(
    t0 = lta_t0,
    t1 = lta_t1
  ) %>%
  map(
    ~ .x %>%
      SummaryTable(
        keepCols = c(
          "Title",
          "LL",
          "Entropy",
          "BIC",        # most trusted stat
          "aBIC",       # remove from kbl
          "T11_KM1LL",
          "T11_VLMR_PValue",
          "BLRT_KM1LL", # remove from kbl
          "BLRT_PValue" # remove from kbl
        ),
        sortBy = "Title"
      ) %>%
      mutate(
        Classes = str_c(str_extract(Title, "\\d"), " Classes"), .before = 1,
        across(where(is.numeric) & !ends_with("Value"), ~ round(.x, 2)),
        across(ends_with("Value"), ~ round(.x, 3)),
        ## number of classes
        ## K = 2:4,
        ## BAYES FACTOR comparision of models with k vs k+1 classes
        BF = signif(exp(-0.5 * BIC - -0.5 * lead(BIC)), 2),
        ## SCHWARTZ INFOMRATION CRITERION (SIC)
        SIC = -0.5 * BIC,
        cmP = exp(SIC - max(SIC))
      ) %>%
      select(-Title, -SIC) %>%
      mutate(
        BF_fx = case_when(
          BF < 3 ~ "weak",
          BF > 3 & BF < 10 ~ "mod",
          BF > 10 ~ "strong",
          TRUE ~ NA_character_
        ),
        cmP_fx = if_else(cmP == max(cmP), "pref", "no_pref")
      ) %>%
      relocate(BF, .before = BF_fx) %>%
      relocate(cmP, .before = cmP_fx) %>%
      filter(!Classes %in% c("1 Classes", "7 Classes"))
  ) %>%
  bind_rows(., .id = "Sample")

.lta_class_select_N <-
  list(
    t0 = lta_t0,
    t1 = lta_t1
  ) %>%
  map(
    `[`, (2)
  ) %>%
  unlist(recursive = FALSE) %>%
  map(`[`, "savedata") %>%
  map(
    ~ .x %>%
      as.data.frame() %>%
      nrow()
  ) %>%
  bind_rows() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(
    Sample = 1,
    N = 2
  ) %>%
  separate(Sample, c("Sample", "NOPA"), sep = "\\.") %>%
  select(-NOPA)

lta_class_select <-
  .lta_class_select %>%
  full_join(.lta_class_select_N) %>%
  relocate(N, .after = Sample)

lta_class_plots <-
  list(
    t0 = lta_t0,
    t1 = lta_t1
  ) %>%
  map(`[`, (2:5)) %>%
  map(
    ~ .x %>%
      plotMixtures(
        coefficients = "stdyx.standardized")
  ) %>%
  wrap_plots()

ggsave(
  "./figures/lta_class_plots.png",
  lta_class_plots
)

### ============================================================================
### END SCRIPT: 1_OPTSEED
### ============================================================================
