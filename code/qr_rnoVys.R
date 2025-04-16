### ============================================================================
### SOURCE SCRIPTS
### ============================================================================
rm(list = ls(all = TRUE))
source("./code/func.R")
source("./code/load.R")
load("./data/proc/lpa_rda.RData")

### ============================================================================
### LOAD DATA
### ============================================================================
## mplus list of candidates with non-zero suicidal ideation at admission
mplus_list <-
  readRDS(here(.m2_lpa_si2))

## final data set with three class lpa profiles added
lpa_fin <-
  readRDS(here(.lpa_fin))

## CREATE DATASETS
## =============================================================================
## LITTLE'S MCAR TEST
## -----------------------------------------------------------------------------
## Test whether assumptions for MI have been met
lpa_mcar <-
  lpa_fin %>%
  filter(session %in% 1 | final_session %in% 1) %>%
  select(
    ID,
    age, sex,
    protocol, session, final_session, sp_db,
    starts_with("bhs_") & ends_with("r"),
    starts_with("bhsSF_") & ends_with("lessone"),
    inq_1:inq_6,
    inq_7r:inq_15r, inq_9, inq_11:inq_12,
    acss_8r:acss_13r, acss_7, acss_11, acss_14, acss_19,
    starts_with(c("dass", "isas_8.")),
    -isas_8.40des,
    mssi_1:mssi_4,
    atsi, cald, hxtrauma, hxsuicide, hxdsh, suicide_attempts, pastpsyc
  ) %>%
  group_by(ID) %>%
  mutate(
    age = min(age, na.rm = TRUE),
    atsi = case_when(
      ## correct atsi status from old sp_db
      sp_db %in% "old" & atsi %in% 4 ~ 0,
      ## collapse categories as cell sizes too small
      atsi %in% 0 ~ 0,
      atsi %in% c(1, 2, 3) ~ 1,
      TRUE ~ NA_real_
    ),
    hxsa = case_when(
      suicide_attempts <= 0          ~ 0,
      suicide_attempts >= 1          ~ 1,
      TRUE                           ~ NA_real_
    ),
    ) %>%
  ungroup() %>%
  select(-suicide_attempts, -sp_db) %>%
  pivot_wider(
    id_cols = c(ID, age, sex, protocol),
    names_from = final_session,
    values_from = c(session, bhs_2r:last_col())
  ) %>%
  select(-ID, -protocol) %>%
  select(!contains("session")) %>%
  naniar::mcar_test() %>%
  rename(chisq = 1, p = 3) %>%
  select(1:3) %>%
  mutate(
    across(
      c(chisq, p),
      ~ sprintf("%.2f", round(.x))
    ),
    df = as.character(df)
  )

## MULTIPLE IMPUTATION
## -----------------------------------------------------------------------------
.lpa_t0t1 <-
  lpa_fin %>%
  filter(!is.na(c)) %>% # filter no class var
  group_by(ID) %>%
  filter(session %in% 1 | final_session %in% 1) %>%
  ungroup() %>%
  mutate(
    session = if_else(is.na(session), -999, session),
    drop = if_else(session %in% 1 & final_session %in% 1, 1, 0), .after = ID,
    across(c(isas_1i, isas_1m), ~ str_extract(.x, "\\d*")),
    across(starts_with("isas_1"), as.numeric),
    isas_intra = isas_intra / 5,
    isas_social = isas_social / 8,
    atsi = case_when(
      ## correct atsi status from old sp_db
      sp_db %in% "old" & atsi %in% 4 ~ 0,
      ## collapse categories as cell sizes too small
      atsi %in% 0 ~ 0,
      atsi %in% c(1, 2, 3) ~ 1,
      TRUE ~ NA_real_
    ),
    hxsa = case_when(
      suicide_attempts <= 0 ~ 0,
      suicide_attempts >= 1 ~ 1,
      TRUE ~ NA_real_
    ),
  ) %>%
  rowwise() %>%
  mutate(dsh_n = sum(c_across(starts_with("isas_1")) > 0, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(drop %in% 0) %>% # select those with more than one session
  ## select vars for mi
  select(
    c, ID, age, sex, session, final_session,
    planned_final, pb, ac, bhs, tb, dep, anx, str, si,
    isas_intra, isas_social, atsi, cald, hxtrauma, hxsuicide, hxdsh, hxsa, pastpsyc, # nolint
    isas_3a, # age onset
    dsh_n
  ) %>%
  group_by(ID) %>%
  mutate(
    age = min(age), planned_final = max(planned_final),
  ) %>%
  ungroup()

.lpa_t0t1_mi <-
  .lpa_t0t1 %>%
  panel_data(id = ID, wave = final_session) %>%
  widen_panel(separator = "_", ignore.attributes = FALSE, varying = NULL) %>%
  ungroup() %>%
  as.data.frame()

lpa_t0t1_mi <-
  .lpa_t0t1_mi %>%
  select(-ac_1, -dsh_n_1) %>%
  rename(ac = ac_0, dsh_n = dsh_n_0) %>%
  ## use age at intake, drop second ac as not routinely collected
  mutate(
    session_1 = if_else(session_1 %in% -999, NA_real_, session_1),
    across(
      c(
        c, ID, sex, session_0, session_1,
        planned_final, atsi, cald, hxtrauma, hxsuicide, hxdsh, hxsa, pastpsyc,
        isas_3a
      ),
      ~ as.integer(as.character(.x))
    )
  )

imp_start <- mice(lpa_t0t1_mi, maxit = 0)

## set method and pred
meth <- imp_start$method
pred <- imp_start$predictorMatrix
pred[, colnames(pred) %in% c("ID")] <- 0 # do not use ID as predictor # nolint
## pred[, colnames(pred) %in% c("ID", "c")] <- 0 # do not use ID or c as predictors

imp_t0t1_100 <-
  parlmice(
    lpa_t0t1_mi,
    maxit = 5,
    meth = meth,
    predictorMatrix = pred,
    cluster.seed = 46810,
    n.core = 4,
    n.imp.core = 25,
    cl.type = "FORK"
  )

lpa_mi_df <-
  imp_t0t1_100 %>%
  complete(., action = "long", include = TRUE) %>%
  group_by(ID) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  mutate(
    si_ch = si_0 - si_1,
    pb_ch = pb_0 - pb_1,
    tb_ch = tb_0 - tb_1,
    bhs_ch = bhs_0 - bhs_1,
    si_ch_r = qnorm((rank(si_ch, na.last = "keep", ties.method = "random") - 0.5) / length(si_ch)), # nolint
    # mice converted this to dbl so back to int
    across(
      c(ID, c, sex, planned_final, atsi:pastpsyc),
      ~ as.factor(as.integer(.x))
    )
  ) %>%
  select(-c(.imp, .id))

lpa_mi_df_long <-
  lpa_mi_df %>%
  long_panel(
    id = "ID",
    wave = "final_session",
    label_location = "end",
    prefix = "_",
    suffix = NULL,
    begin = 0,
    end = 1
  ) %>%
  ungroup() %>%
  mutate(final_session = as.factor(final_session))

lpa_t0t1_outc_w <-
  lpa_mi_df

## SUBSET TIME 1: DEMO
## -----------------------------------------------------------------------------
lpa_t0_og <- # orginal dataset, no MI
  lpa_fin %>%
  filter(!is.na(c), if_any(c(session, final_session), ~.x %in% 1)) %>%
  group_by(ID) %>%
  mutate(
    total_sessions = max(session), .after = session,
    planned_final = max(planned_final)
  ) %>%
  ungroup() %>%
  filter(session %in% 1) %>%
  select(ID:suburb)

lpa_t0 <-
  lpa_mi_df_long %>%
  filter(!is.na(c), if_any(c(session, final_session), ~.x %in% 1)) %>%
  group_by(ID) %>%
  mutate(
    total_sessions = max(session), .after = session,
    planned_final = max(as.numeric(planned_final))
  ) %>%
  ungroup() %>%
  filter(session %in% 1)

## SUBSET TIME 1 WITH MATCHING TIME 2
## -----------------------------------------------------------------------------
## TODO remove unused code below
lpa_t0t1 <-
  lpa_mi_df_long %>%
  filter(!is.na(c), if_any(c(session, final_session), ~.x %in% 1)) %>%
  select(
    where(is.numeric) & !contains("ch"),
    c, ID, sex, final_session, planned_final
  )

## LPA WITH FINAL SESSION BY CLASS
## -----------------------------------------------------------------------------
## here we look at t0 to t1 change
.lpa_t0t1_outc <-
  lpa_mi_df_long %>%
  filter(!is.na(c), if_any(c(session, final_session), ~.x %in% 1)) %>%
  mutate(sex = factor(if_else(sex %in% 2, NA_character_, as.character(sex)))) %>% # nolint
  group_by(ID) %>%
  mutate(total_sessions = max(session), .after = session) %>%
  mutate(
    cov_sess = max(session), .after = final_session,
    ac = ac[session %in% 1]
  ) %>%
  ungroup()

## create difference df for change score reliability
.lpa_alt_dif <-
  lpa_fin %>%
  filter(session %in% 1 | final_session %in% 1) %>%
  group_by(ID) %>%
  mutate(across(starts_with(c("isas_", "is", "dsh_n", "diag")), ~ max(.x))) %>%
  ungroup() %>%
  panel_data(id = ID, wave = final_session) %>%
  widen_panel(separator = "_", ignore.attributes = FALSE, varying = NULL) %>%
  ungroup() %>%
  as.data.frame()

lpa_alt_dif <-
  .lpa_alt_dif %>%
  mutate(
    across(c(inq_1_0:inq_6_0), .names = "{.col}_ch") -
      across(c(inq_1_1:inq_6_1)),
    across(c(inq_7r_0:inq_15r_0, inq_9_0, inq_11_0:inq_12_0), .names = "{.col}_ch") - # nolint
      across(c(inq_7r_1:inq_15r_1, inq_9_1, inq_11_1:inq_12_1)),
    across(c(mssi_1_0:mssi_4_0, starts_with("rc_mssi") & ends_with("_0")), .names = "{.col}_ch") - # nolint
      across(c(mssi_1_1:mssi_4_1, starts_with("rc_mssi") & ends_with("_1"))),
    ## get bhs sum of protocol 1
    across(c(starts_with("bhs_") & ends_with("r_0")), .names = "{.col}_ch") -
      across(c(starts_with("bhs_") & ends_with("r_1"))),
    ## get bhs sum of protocol 2
    across(c(starts_with("bhsSF_") & ends_with("lessone_0")), .names = "{.col}_ch") - # nolint
      across(c(starts_with("bhsSF_") & ends_with("lessone_1"))) # nolint
  ) %>%
  rename_with(~ str_remove(.x, "_0"), ends_with("_0_ch")) %>%
  select(
    c, ID, sex, starts_with(c("session_", "age_")),
    protocol,
    inq_1_0:inq_6_0, inq_1_1:inq_6_1,
    inq_7r_0:inq_15r_0, inq_9_0, inq_11_0:inq_12_0,
    inq_7r_1:inq_15r_1, inq_9_1, inq_11_1:inq_12_1,
    mssi_1_0:mssi_4_0, starts_with("rc_mssi") & ends_with("_0"),
    mssi_1_1:mssi_4_1, starts_with("rc_mssi") & ends_with("_1"),
    starts_with("bhs_") & ends_with("r_0"),
    starts_with("bhs_") & ends_with("r_1"),
    starts_with("bhsSF_") & ends_with("lessone_0"),
    starts_with("bhsSF_") & ends_with("lessone_1"),
    ends_with("ch"),
  ) %>%
  .pro_mean("bhs_m1_ch", starts_with("bhs_") & ends_with("r_ch")) %>%
  ## get bhs sum of protocol 2
  .pro_mean("bhs_m2_ch", starts_with("bhsSF_") & ends_with("lessone_ch")) %>%
  .pro_mean("pb_ch", inq_1_ch:inq_6_ch) %>%
  .pro_mean("tb_ch", c(inq_7r_ch:inq_15r_ch, inq_9_ch, inq_11_ch:inq_12_ch)) %>%
  ## caolesce protocol 1 and 2 into one score
  mutate(bhs_ch = coalesce(bhs_m1_ch, bhs_m2_ch)) %>%
  .pro_mean("si_ch", c(mssi_1_ch:mssi_4_ch, starts_with("rc_mssi")) & ends_with("_ch")) %>% # nolint
  long_panel(
    prefix = "_",
    periods = c("0", "1", "ch"),
    id = "ID",
    wave = "final_session",
    label_location = "end"
  ) %>%
  ungroup()

lpa_t0t1_outc <-
  .lpa_t0t1_outc %>%
  select(
    c, ID, age, sex, session, final_session, cov_sess, planned_final, pb, ac,
    bhs, tb, dep, anx, str, si, isas_intra, isas_social, dsh_n,
    hxsa, hxsuicide, hxdsh, hxtrauma, atsi, cald,
  ) %>%
  group_by(final_session) %>%
  mutate(
    across(
      c(
        si, pb, tb, bhs, dep, anx, str, ac, isas_social, isas_intra,
        age,
      ),
      ~ scale(.x),
      .names = "{.col}_z"
    ),
  ) %>%
  ungroup() %>%
  mutate(
    fin_ses_c = with(., factor(final_session):factor(c)), .after = c,
    across(c(c, sex, ID, final_session, planned_final), ~ as.factor(.x))
  )

## SECONDARY ANALYSIS BY C DATA
## -----------------------------------------------------------------------------
## recode values of demo vars to make fu df
.fu <-
  .lpa_t0t1_outc %>% # TODO use MI data set
  filter(session %in% 1) %>%
  mutate(
    sex = case_when(
      sex %in% 0 ~ "Female",
      sex %in% 1 ~ "Male",
      sex %in% 2 ~ "Other",
      TRUE ~ NA_character_
    ),
    across(
      c(planned_final, atsi:pastpsyc),
      ~ factor(case_when(
        .x %in% 0 ~ "No",
        .x %in% 1 ~ "Yes",
        TRUE ~ as.character(.x)
      )),
      .names = "{.col}2"
    ),
    across(c(c, sex), as.factor)
  )
.fu

### ============================================================================
### ANALYSESES
### ============================================================================

## LPA CLASS TABL
## =============================================================================
lpa_class_tab <-
  SummaryTable(
    mplus_list,
    keepCols = c(
      "Title",
      "LL",
      "Entropy",
      ## "AIC",
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
    across(ends_with("Value"), ~ round(.x, 3))
  ) %>%
  select(-Title) %>%
  filter(!Classes %in% "1 Classes")
lpa_class_tab

## RELIABILITY
## =============================================================================
.pb_alpha <-
  lpa_alt_dif %>%
  filter(!is.na(c)) %>%
  select(final_session, inq_1:inq_6)
.tb_alpha <-
  lpa_alt_dif %>%
  filter(!is.na(c)) %>%
  select(final_session, inq_7r:inq_15r, inq_9, inq_11:inq_12)
.bhs_20 <-
  lpa_alt_dif %>%
  filter(!is.na(c)) %>%
  filter(protocol %in% 1) %>%
  select(final_session, starts_with("bhs_") & ends_with("r"))
.bhs_sf <-
  lpa_alt_dif %>%
  filter(!is.na(c)) %>%
  filter(protocol %in% 2) %>%
  select(final_session, starts_with("bhsSF_") & ends_with("lessone"))
.si_alpha <-
  lpa_alt_dif %>%
  filter(!is.na(c)) %>%
  select(final_session, mssi_1:mssi_4, starts_with("rc_mssi"))
.si_alpha_4 <-
  lpa_alt_dif %>%
  filter(!is.na(c)) %>%
  select(final_session, mssi_1:mssi_4)
.ac_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, c(acss_8r:acss_13r, acss_7, acss_11, acss_14, acss_19))
.dep_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1 | final_session %in% 1) %>%
  select(final_session, num_range("dass_", c(3, 5, 10, 13, 16, 17, 21)))
.anx_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1 | final_session %in% 1) %>%
  select(final_session, num_range("dass_", c(2, 4, 7, 9, 15, 19, 20)))
.str_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1 | final_session %in% 1) %>%
  select(final_session, num_range("dass_", c(1, 6, 8, 11, 12, 14, 18)))
.is_aff_reg_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(1, 14, 27)))
.is_int_bnd_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(2, 15, 28)))
.is_slf_pun_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(3, 16, 29)))
.is_slf_car_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(4, 17, 30)))
.is_ant_dis_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(5, 18, 31)))
.is_ant_sui_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(6, 19, 32)))
.is_sen_see_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(7, 20, 33)))
.is_peer_bnd_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(8, 21, 34)))
.is_int_inf_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(9, 22, 35)))
.is_tough_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(10, 23, 36)))
.is_mrk_dis_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(11, 24, 37)))
.is_revenge_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(12, 25, 38)))
.is_autnmy_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(13, 36, 39)))
.isas_intra_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(
    1, 3, 5, 6, 11, 14, 16, 18,
    19, 24, 27, 29, 31, 32, 37
  )))
.isas_social_alpha <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  filter(session %in% 1) %>%
  select(final_session, num_range("isas_8.", c(
      2, 4, 7, 8, 9, 10, 12, 13, 15, 17,
      20, 21, 22, 23, 25, 26, 28, 30,
      33, 34, 35, 36, 38, 39
  )))

## create list of reliability stats for measures
meas_rel <-
  list(
    pb = .pb_alpha,
    tb = .tb_alpha,
    bhs_20 = .bhs_20,
    bhs_sf = .bhs_sf,
    si = .si_alpha,
    si_4 = .si_alpha_4,
    ac = .ac_alpha,
    dep = .dep_alpha,
    anx = .anx_alpha,
    str = .str_alpha,
    is_aff_reg = .is_aff_reg_alpha,
    is_int_bnd = .is_int_bnd_alpha,
    is_slf_pun = .is_slf_pun_alpha,
    is_slf_car = .is_slf_car_alpha,
    is_ant_dis = .is_ant_dis_alpha,
    is_ant_sui = .is_ant_sui_alpha,
    is_sen_see = .is_sen_see_alpha,
    is_peer_bnd = .is_peer_bnd_alpha,
    is_int_inf = .is_int_inf_alpha,
    is_tough = .is_tough_alpha,
    is_mrk_dis = .is_mrk_dis_alpha,
    is_revenge = .is_revenge_alpha,
    is_autnmy = .is_autnmy_alpha,
    isas_intra = .isas_intra_alpha,
    isas_social = .isas_social_alpha
  ) %>%
  map(
    ~ .x %>%
      split(.$final_session) %>%
      map(
        ~ .x %>%
          select(-final_session) %>%
          as.data.frame() %>%
          psych::alpha()
      )
  ) %>%
  unlist(recursive = FALSE) %>%
  map(`[`, "total") %>%
  unlist(recursive = FALSE) %>%
  map(`[`, "std.alpha") %>%
  unlist() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(val = 2) %>%
  mutate(
    rowname = str_remove(rowname, "\\.total\\.std\\.alpha"),
  ) %>%
  separate(rowname, c("var", "final_session"), sep = "\\.") %>%
  pivot_wider(
    names_from = final_session,
    values_from = val,
    names_prefix = "fin_"
  ) %>%
  mutate(
    across(
      starts_with("fin_"),
      ~ str_remove(sprintf("%.2f", round(.x, 2)), "^0")
    )
  )

## ANOVAS BY C
## =============================================================================
## TODO get the standardised residuals of each model
## https://www.statology.org/standardized-residuals-in-r/
## NOTE the below is a bit superfluous as it doubles the follow up analysis
.lpa_anova <-
  lpa_t0t1_outc_w %>%
  select(ID, c, !ends_with(c("_z", "_r")) & where(is.numeric), -session_0) %>%
  rename(session_0 = session_1) %>%
  pivot_longer(
    cols = !c(ID:c),
    names_to = c("name", "time"),
    names_sep = "_"
  ) %>%
  mutate(
    name = if_else(name %in% c("isas", "dsh"), str_c(name, "_", time), name),
    time = if_else(!time %in% c(0:1, "ch"), "0", time)
  ) %>%
  split(.$time)


lpa_md_ci <-
  list(
    `2v1` =
      .lpa_anova %>%
        map(
          ~ .x %>%
            mutate(c = as.character(c)) %>%
            filter(!c %in% 3) %>%
            pivot_wider(
              names_from = c, values_from = value,
              names_glue = "c{.name}"
            ) %>%
            split(.$name) %>%
            map(
              ~ .x %$%
                DescTools::MeanDiffCI(c2, c1, na.rm = TRUE)
            ) %>%
            bind_rows(.id = "var")
        ) %>%
        set_names(c("t0", "t1", "ch")),
    `3v1` =
      .lpa_anova %>%
        map(
          ~ .x %>%
            mutate(c = as.character(c)) %>%
            filter(!c %in% 2) %>%
            arrange(desc(c)) %>%
            pivot_wider(
              names_from = c, values_from = value,
              names_glue = "c{.name}"
            ) %>%
            split(.$name) %>%
            map(
              ~ .x %$%
                DescTools::MeanDiffCI(c3, c1, na.rm = TRUE)
            ) %>%
            bind_rows(.id = "var")
        ) %>%
        set_names(c("t0", "t1", "ch")),
    `2v3` =
      .lpa_anova %>%
        map(
          ~ .x %>%
            mutate(c = as.character(c)) %>%
            filter(!c %in% 1) %>%
            arrange(c) %>%
            pivot_wider(
              names_from = c, values_from = value,
              names_glue = "c{.name}"
            ) %>%
            split(.$name) %>%
            map(
              ~ .x %$%
                DescTools::MeanDiffCI(c2, c3, na.rm = TRUE)
            ) %>%
            bind_rows(.id = "var")
        ) %>%
        set_names(c("t0", "t1", "ch"))
  ) %>%
  map(~ .x %>% bind_rows(.id = "time")) %>%
  bind_rows(.id = "cont") %>%
  rename(
    md = meandiff,
    ci_l = lwr.ci,
    ci_u = upr.ci
  ) %>%
  pivot_wider(
    names_from = cont,
    values_from = md:ci_u,
  )

lpa_d_ci <-
  .lpa_anova %>%
  map(
    ~ .x %>%
      split(.$name) %>%
      map(
        ~ .x %>%
          group_by(c) %>%
          summarise(
            n = n(),
            across(
              value,
              list(
                m = ~ mean(.x, na.rm = TRUE),
                sd = ~ sd(.x, na.rm = TRUE)
              )
            )
          ) %>%
          rename_with(~ str_remove(.x, "value_")) %>%
          pivot_wider(
            names_from = c,
            values_from = n:last_col()
          )
      ) %>%
      bind_rows(.id = "var") %>%
      mutate(
        d_2v1 = compute.es::mes(m_2, m_1, sd_2, sd_1, n_2, n_1, verbose = FALSE)$d, # nolint
        d_2v1_l = compute.es::mes(m_2, m_1, sd_2, sd_1, n_2, n_1, verbose = FALSE)$l.d, # nolint
        d_2v1_u = compute.es::mes(m_2, m_1, sd_2, sd_1, n_2, n_1, verbose = FALSE)$u.d, # nolint
        d_3v1 = compute.es::mes(m_3, m_1, sd_3, sd_1, n_3, n_1, verbose = FALSE)$d, # nolint
        d_3v1_l = compute.es::mes(m_3, m_1, sd_3, sd_1, n_3, n_1, verbose = FALSE)$l.d, # nolint
        d_3v1_u = compute.es::mes(m_3, m_1, sd_3, sd_1, n_3, n_1, verbose = FALSE)$u.d, # nolint
        d_2v3 = compute.es::mes(m_2, m_3, sd_2, sd_3, n_2, n_3, verbose = FALSE)$d, # nolint
        d_2v3_l = compute.es::mes(m_2, m_3, sd_2, sd_3, n_2, n_3, verbose = FALSE)$l.d, # nolint
        d_2v3_u = compute.es::mes(m_2, m_3, sd_2, sd_3, n_2, n_3, verbose = FALSE)$u.d, # nolint
      )
  ) %>%
  set_names(c("t0", "t1", "ch")) %>%
  bind_rows(.id = "time")

lpa_an <-
  .lpa_anova %>%
  map(
    ~ .x %>%
      split(.$name) %>%
      map(~ tidy(Anova(lm(value ~ c, .x), type = 3))) %>%
      bind_rows(.id = "var") %>%
      filter(term %in% "c") %>%
      select(var, statistic, p.value)
  ) %>%
  set_names(c("t0", "t1", "ch")) %>%
  bind_rows(.id = "time")

lpa_tt <-
 .lpa_anova %>%
  map(
    ~ .x %>%
      split(.$name) %>%
      map(
        ~ .x %$%
          pairwise.t.test.with.t.and.df(value, c, p.adjust.method = "bonf")
      ) %>%
      map(`[`, c("p.value", "t.value", "dfs")) %>%
      map(~ as.data.frame(.x) %>%
        clean_names() %>%
        rownames_to_column("cont")) %>%
      bind_rows(.id = "var")
  ) %>%
  bind_rows(.id = "time") %>%
  pivot_longer(
    cols = !c(time, var, cont, dfs),
    names_to = c("test", "thing", "cont2"),
    names_pattern = "^(p|t)(.value.)([1-2])",
    values_drop_na = TRUE
  ) %>%
  select(-thing) %>%
  pivot_wider(names_from = test, values_from = value) %>%
  relocate(dfs, .after = last_col()) %>%
  split(.$time) %>%
  set_names(c("t0", "t1", "ch")) %>%
  map(
    ~ .x %>%
      unite("contr", c(cont, cont2), sep = "v") %>%
      pivot_wider(names_from = contr, values_from = c(t, p))
  ) %>%
  bind_rows(.id = "time")

lpa_anova_fx_mes <-
  lpa_an %>%
  full_join(lpa_tt) %>%
  full_join(lpa_md_ci) %>%
  full_join(lpa_d_ci) %>%
  mutate(
    across(
      starts_with("d_"),
      ~ case_when(
        abs(.x) >= 0 & abs(.x) < 0.5 ~ "small",
        abs(.x) >= 0.5 & abs(.x) < 0.8 ~ "medium",
        abs(.x) >= 0.8 ~ "large"
      ),
      .names = "{.col}_fx"
    ),
    var = case_when(
      var %in% "isas_intra" ~ "intra",
      var %in% "isas_social" ~ "social",
      var %in% "dsh_n" ~ "DSHn",
      var %in% "isas_3a" ~ "DSHage",
      TRUE ~ var
    )
  ) %>%
  split(.$time)

lpa_sample <-
  lpa_t0 %>%
  select(c, pb, tb, bhs, ac) %>%
  mutate(
    across(
      c(pb, tb, bhs, ac),
      ~ scale(.x),
      .names = "{.col}_z"
    ),
  ) %>%
  group_by(c) %>%
  summarise(
    N = n(),
    across(
      c(pb, tb, bhs, ac),
      list(
        m = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      )
    ),
    across(
      ends_with(c("_z")),
      ~ round(mean(.x, na.rm = TRUE), 2)
    )
  ) %>%
  mutate(
    prc = N / sum(N) * 100, .after = N
  ) %>%
  select(c, N, prc, ends_with(c("_m", "_sd", "_z")))

## FOLLOW UP ANALYSIS BY C
## -----------------------------------------------------------------------------
## CONTINUOUS VARS
## -------------------------------
## CLASS SIZES FULL ONLY
.lpa_n_fu_full <-
  .fu %>%
  select(ID) %>%
  summarise(n_sd = sum(!is.na(ID))) %>%
  mutate(n_m = round(n_sd / sum(n_sd) * 100, 2), c = "full")

## CLASS SIZES ALL
.lpa_n_fu <-
  .fu %>%
  select(c, ID) %>%
  group_by(c) %>%
  summarise(n_sd = sum(!is.na(ID))) %>%
  mutate(
    n_m = round(n_sd / sum(n_sd) * 100, 2),
    c = as.character(c)
  ) %>%
  bind_rows(.lpa_n_fu_full) %>%
  pivot_longer(
    cols = !c,
    names_to = "key",
    values_to = "value"
  ) %>%
  mutate(
    var = str_replace_all(key, "_m|_sd", ""),
    type = str_extract(key, "(m|sd)$")
  ) %>%
  pivot_wider(
    names_from = c(type, c),
    values_from = value,
  ) %>%
  select(-key) %>%
  group_by(var) %>%
  summarise(across(everything(), .coalesce_by_column)) %>%
  ungroup() %>%
  mutate(var = str_replace(var, "n", "Class size")) %>%
  select(var, starts_with(c("m", "sd"))) %>%
  relocate(ends_with("_full"), .after = var)

## UNIVARIATE ANOVAS OF COVARIATES
.fu_test <-
  .fu %>%
  select(
    pb, tb, bhs, ac,
    si, dep, anx, str, age,
    isas_intra, isas_social, total_sessions,
    dsh_n, isas_3a,
  ) %>%
  map(~ tidy(Anova(lm(. ~ .fu$c, na.action = na.exclude), type = 3))) %>%
  unlist(recursive = FALSE) %>%
  as_tibble() %>%
  filter(if_any(ends_with("term"), ~ .x %in% ".fu$c")) %>%
  select(ends_with(c("statistic", "p.value"))) %>%
  pivot_longer(
    everything(),
    names_to = c("var", "stat"),
    names_sep = "\\."
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  )

## INDIVIDUAL CONTRASTS OF FU COVARIATE DATA
.fu_test_cntrst <-
  .fu %>%
  select(
    pb, tb, bhs, ac,
    si, dep, anx, str, age,
    isas_intra, isas_social, total_sessions,
    dsh_n, isas_3a,
  ) %>%
  mutate(across(where(is.character), ~ as.numeric(.x))) %>%
  map(~ tidy(pairwise.t.test(., .fu$c, p.adjust.method = "bonf"))) %>%
  unlist(recursive = FALSE) %>%
  as_tibble() %>%
  mutate(across(where(is.character), as.numeric)) %>%
  pivot_longer(
    col = contains("value"),
    names_to = "key",
    values_to = "value"
  ) %>%
  select(contains(c("group1", "group2")), everything()) %>%
  arrange(key) %>%
  ## TODO find tidy way to do this
  .[!duplicated(lapply(., summary))] %>%
  rename(group1 = 1, group2 = 2, var = 3, p_value = 4) %>%
  mutate(var = str_replace(var, "\\.p\\.value", "")) %>%
  pivot_wider(
    names_from = c(group1, group2),
    values_from = p_value
  ) %>%
  clean_names()

## CONTINOUS MEANS
.lpa_m_sd_fu_full <-
  .fu %>%
  select(
    pb, tb, bhs, ac,
    si, dep, anx, str, age,
    isas_intra, isas_social, total_sessions,
    dsh_n, isas_3a,
  ) %>%
  mutate(across(everything(), as.numeric)) %>%
  summarise(
    across(
      everything(),
      list(
        m = ~ round(mean(.x, na.rm = TRUE), 2),
        sd = ~ round(sd(.x, na.rm = TRUE), 2)
      )
    )
  ) %>%
  mutate(c = "full", .before = 1) %>%
  pivot_longer(
    cols = !c,
    names_to = "key",
    values_to = "value"
  ) %>%
  mutate(
    var = str_replace_all(key, "(?!_mrk_)(_m|_sd)", ""),
    type = str_extract(key, "(m|sd)$"),
    ) %>%
  pivot_wider(
    names_from = c(type, c),
    values_from = value,
    ) %>%
  select(-key) %>%
  group_by(var) %>%
  summarise(across(everything(), .coalesce_by_column)) %>%
  ungroup()

## MEANS AND SD BY CLASS
.lpa_m_sd_fu <-
  .fu %>%
  select(
    c,
    pb, tb, bhs, ac,
    si, dep, anx, str, age,
    isas_intra, isas_social, total_sessions,
    dsh_n, isas_3a,
  ) %>%
  mutate(across(everything(), as.numeric)) %>%
  group_by(c) %>%
  summarise(
    across(
      everything(),
      list(
        m = ~ round(mean(.x, na.rm = TRUE), 2),
        sd = ~ round(sd(.x, na.rm = TRUE), 2)
      )
    )
  ) %>%
  ungroup() %>%
  pivot_longer(
    cols = !c,
    names_to = "key",
    values_to = "value"
  ) %>%
  mutate(
    var = str_replace_all(key, "(?!_mrk_)(_m|_sd)", ""),
    type = str_extract(key, "(m|sd)$"),
    ) %>%
  pivot_wider(
    names_from = c(type, c),
    values_from = value,
    ) %>%
  select(-key) %>%
  group_by(var) %>%
  summarise(across(everything(), .coalesce_by_column)) %>%
  ungroup() %>%
  full_join(.lpa_m_sd_fu_full) %>%
  filter(!var %in% "c") %>%
  select(var, ends_with("full"), everything())

## FULL CONTINUOUS FU TABLE
fu_test_cont <-
  .lpa_n_fu %>%
  full_join(.lpa_m_sd_fu) %>%
  full_join(.fu_test) %>%
  full_join(.fu_test_cntrst) %>%
  mutate(test = "cont", bonf = 0.05)

fu_test_cont %>%
  select(!ends_with("full")) %>%
  filter(p < .05) %>%
  filter(if_any(c(x2_1:x3_2), ~.x < .05)) %>%
  print(n = Inf)

## CATEGORICAL FU VARS
## --------------------------
## CHISQ ON FU COVARIATE DATA
.fu_test_nonp <-
  .fu %>%
  select(sex, ends_with("2")) %>%
  map(~ tidy(chisq.test(., .fu$c))) %>%
  unlist(recursive = FALSE) %>%
  as_tibble() %>%
  filter(if_any(ends_with("term"), ~ .x %in% ".fu$c")) %>%
  select(ends_with(c("statistic", "p.value"))) %>%
  pivot_longer(everything(), names_to = c("var", "stat"), names_sep = "\\.") %>%
  pivot_wider(names_from = stat, values_from = value)

.fu_test_nonp %>% print(n = Inf)

## GET STD RESIDUALS
.fu_resd_use <-
  .fu %>%
  select(sex, ends_with("2")) %>%
  map(~ chisq.test(table(., .fu$c))) %>%
  map(`[`, c("stdres", "statistic", "p.value")) %>%
  map(~ as.data.frame(.)) %>%
  map(~ pivot_wider(., names_from = 2, values_from = 3)) %>%
  map(~ janitor::clean_names(.)) %>%
  map(.,
      ~ as.data.frame(bind_cols(
        .x,
        ## create bonf corrected cutoff
        bonf = qnorm(1 - 0.05 / (3 * nrow(.x)) / 2),
        ))
      ) %>%
  bind_rows(., .id = "var") %>%
  rename(
    ## to match cont df
    value = stdres,
    p = p_value,
    x2_1 = x1,
    x3_1 = x2,
    x3_2 = x3
  )
.fu_resd_use

## CREATE POSTHOC MULTIPLE COMPARISION GROUPS
## 2v1
.chiq_bonf_21 <-
  .fu %>%
  select(sex, ends_with("2"), c) %>%
  mutate(
    across(where(is.factor), ~ as.character(.x))
  ) %>%
  filter(!c %in% 3)
## 3v1
.chiq_bonf_31 <-
  .fu %>%
  select(sex, ends_with("2"), c) %>%
  mutate(
    across(where(is.factor), ~ as.character(.x))
  ) %>%
  filter(!c %in% 2)
## 2v3
.chiq_bonf_23 <-
  .fu %>%
  select(sex, ends_with("2"), c) %>%
  mutate(
    across(where(is.factor), ~ as.character(.x))
  ) %>%
  filter(!c %in% 1)

## GET MULTIPLE CHISQ COMPARISONS
.fu_chisq_bonf <-
  list(
    x2_1 = .chiq_bonf_21,
    x3_1 = .chiq_bonf_31,
    x3_2 = .chiq_bonf_23
  ) %>%
  map(
    ~ .x %>%
      mutate(across(where(is.character), ~ as.factor(.x))) %>%
      pivot_longer(-c) %>%
      split(.$name) %>%
      map(~ .x %$% tidy(chisq.test(value, c))) %>%
      unlist(recursive = FALSE) %>%
      as_tibble() %>%
      ## select(ends_with(c("statistic", "p.value"))) %>%
      select(ends_with("p.value")) %>%
      pivot_longer(everything(), names_to = c("var", "stat"), names_sep = "\\.") %>% # nolint
      pivot_wider(names_from = stat, values_from = value)
  ) %>%
  bind_rows(.id = "cont") %>%
  pivot_wider(names_from = cont, values_from = p) %>%
  mutate(bonf = 0.05 / 3)

.fu_nonp_use <-
  full_join(.fu_test_nonp, .fu_chisq_bonf)

## GET FULL SAMPLE DATA
.lpa_m_sd_fu_nonp_full <-
  .fu %>%
  select(c, sex, ends_with("2")) %>%
  pivot_longer(cols = -c) %>%
  count(name, value) %>%
  filter(!is.na(value)) %>%
  mutate(c = "full", .before = 1) %>%
  group_by(c, name) %>%
  rename(sd = n) %>%
  mutate(across(
    sd,
    ~ round(.x / sum(.x) * 100, 2),
    .names = "m"
  )) %>%
  mutate(
    row = row_number(),
    across(where(is.integer), as.double)
  ) %>%
  pivot_wider(
    names_from = c,
    values_from = c(sd, m),
  ) %>%
  filter(!name %in% "c") %>%
  select(-row)
.lpa_m_sd_fu_nonp_full %>% print(n = Inf)

fu_test_cont_nonp <-
  .fu %>%
  select(c, sex, ends_with("2")) %>%
  pivot_longer(cols = -c) %>%
  count(c, name, value) %>%
  filter(!is.na(value)) %>%
  group_by(c, name) %>%
  rename(sd = n) %>%
  pivot_wider(names_from = c, names_prefix = "sd_", values_from = sd) %>%
  ## create percentages
  group_by(name) %>%
  mutate(
    across(
      sd_1:sd_3,
      ~ round(.x / sum(.x, na.rm = TRUE) * 100, 2),
      .names = "m_{1:3}"
    )) %>%
  ungroup() %>%
  full_join(.lpa_m_sd_fu_nonp_full) %>%
  arrange(name, value) %>%
  group_by(name, value) %>%
  arrange(name, desc(m_full)) %>%
  ungroup() %>%
  select(name, value, ends_with(c("full", "1", "2", "3")), everything()) %>%
  rename(var = name) %>%
  full_join(.fu_nonp_use) %>%
  ## full_join(.fu_resd_use) %>%
  mutate(test = "fact", .before = bonf)

fu_test_cont_nonp %>% print(n = Inf)
fu_test_cont %>% print(n = Inf)

.fu_test_full <-
  fu_test_cont %>%
  full_join(fu_test_cont_nonp) %>%
  relocate(c(var, value), .before = 1) %>%
  filter(!var %in% "planned_final2") %>%
  mutate(
    order = case_when(
      var %in% "Class size"     ~ 1,
      var %in% "tb"             ~ 2,
      var %in% "pb"             ~ 3,
      var %in% "bhs"            ~ 4,
      var %in% "ac"             ~ 5,
      var %in% "sex"            ~ 6,
      var %in% "age"            ~ 7,
      var %in% "atsi2"          ~ 8,
      var %in% "cald2"          ~ 9,
      var %in% "dep"            ~ 10,
      var %in% "anx"            ~ 11,
      var %in% "str"            ~ 12,
      var %in% "isas_social"    ~ 13,
      var %in% "isas_intra"     ~ 14,
      var %in% "dsh_n"          ~ 15,
      var %in% "isas_3a"        ~ 16,
      var %in% "hxdsh2"         ~ 17,
      var %in% "si"             ~ 18,
      var %in% "hxsa2"          ~ 19,
      var %in% "hxsuicide2"     ~ 20,
      var %in% "hxtrauma2"      ~ 21,
      var %in% "pastpsyc2"      ~ 22,
      var %in% "total_sessions" ~ 23
    ), .before = 1,
    type = case_when(
      order %in% 1     ~ "Class size<sup>a</sup> <i>[%(n)]</i>",
      order %in% 2:5   ~ "LPA indicator",
      order %in% 6:9   ~ "Demographic",
      order %in% 10:12 ~ "Pyschological distress",
      order %in% 13:17 ~ "Deliberate self-injury",
      order %in% 18:21 ~ "Suicide-related",
      TRUE             ~ "Treatment-related"
    ),
    .var2 = var %>% str_remove("\\d$"),
    var3 = var %>% str_remove("\\d$"),
    .value = tolower(as.character(value)),
    value3 = tolower(as.character(value)),
    value3 = if_else(value3 %in% "female", "yfemale", value3),
    var = str_remove_all(str_replace(var, "_", " "), "\\d$"),
    var = case_when(
      str_detect(var, "^(pb|tb|bhs|ac|atsi|cald|dep|anx|str)") ~ toupper(var), # nolint
      str_detect(var, "^age|^total") ~ str_to_sentence(var),
      str_detect(var, "sex") ~ str_to_title(paste0(var, " (", value, ")")),
      var %in% "si" ~ "MSSI",
      str_detect(var, "^Class") ~ NA_character_,
      str_detect(var, "^isas|^hx|^sa|^trauma|^suicide") ~ str_replace_all(
        var, c(
          "^(isas)( social)"     = "ISAS<sub>Inter</sub>",
          "^(isas)( intra)"      = "ISAS<sub>Intra</sub>",
          "^(isas)( 3a)"         = "Age of onset",
          "^(hx)(trauma)(\\<*)"  = "Trauma<sub>Hx</sub> <i>[%(n)]</i>",
          "^(hx)(suicide)(\\<*)" = "SI<sub>Hx</sub> <i>[%(n)]</i>",
          "^(hx)(sa)(\\<*)"      = "SA<sub>Hx</sub> <i>[%(n)]</i>",
          "^(hx)(dsh)(\\<*)"     = "DSH<sub>Hx</sub> <i>[%(n)]</i>"
        )
      ),
      str_detect(var, "dsh n") ~ "Number of methods",
      str_detect(var, "^past") ~ "Previous intervention  <i>[%(n)]</i>",
      str_detect(var, "session ") ~ str_replace(var, "session ", "sessions "),
      TRUE ~ str_to_sentence(var)
    ),
    across(
      starts_with("sd_"),
      ~ case_when(
        .data$test %in% "cont" && !.data$var %in% "Class size" ~
          sprintf("%.2f", round(.x, 2)),
        TRUE ~ as.character(.x)
      )
    ),
    across(starts_with("m_"), ~ sprintf("%.2f", round(.x, 2))),
    full    = str_c(m_full, "(", sd_full, ")"),
    low_si  = str_c(m_1, "(", sd_1, ")"),
    hi_si   = str_c(m_2, "(", sd_2, ")"),
    mod_si  = str_c(m_3, "(", sd_3, ")"),
    chisq_f = case_when(
      ## html tags for supescript
      p < 0.001 ~ paste0(sprintf("%.2f", statistic), "<sup>***</sup>"),
      p < 0.01  ~ paste0(sprintf("%.2f", statistic), "<sup>** </sup>"),
      p < 0.05  ~ paste0(sprintf("%.2f", statistic), "<sup>*  </sup>"),
      p > 0.05  ~ paste0(sprintf("%.2f", statistic), "<sup>   </sup>"),
      TRUE ~ NA_character_
    ),
  ) %>%
  arrange(order, value) %>%
  relocate(c(order, type, var, value), .before = 1) %>%
  unite(var2, c(.var2, .value)) %>%
  mutate(var2 = str_remove(var2, "_NA$"))

## get class numbers
.n_c123 <-
  .fu_test_full %>%
  filter(if_all(var2, ~ str_detect(.x, "^Class"))) %>%
  select(starts_with("sd_") & -sd_full) %>%
  mutate(across(everything(), ~ as.numeric(.x)))

fu_test <-
  .fu_test_full %>%
  split(.$var3) %>%
  map(
    ~ .x %>%
      arrange(value3) %>%
      mutate(
        across(starts_with(c("m_", "sd_")), ~ as.numeric(.x)),
        diff_2v1 = case_when(
          ## odds ratio
          p < 0.05 & test %in% "fact" ~ (sd_2[2] / sd_2[1]) / (sd_1[2] / sd_1[1]),                                           # nolint
          ## mean diff
          p < 0.05 & test %in% "cont" ~ m_2 - m_1,
          TRUE ~ NA_real_
        ),
        diff_3v1 = case_when(
          p < 0.05 & test %in% "fact" ~ (sd_3[2] / sd_3[1]) / (sd_1[2] / sd_1[1]),                                           # nolint
          p < 0.05 & test %in% "cont" ~ m_3 - m_1,
          TRUE ~ NA_real_
        ),
        diff_2v3 = case_when(
          p < 0.05 & test %in% "fact" ~ (sd_2[2] / sd_2[1]) / (sd_3[2] / sd_3[1]),                                           # nolint
          p < 0.05 & test %in% "cont" ~ m_2 - m_3,
          TRUE ~ NA_real_
        )
      )
  ) %>%
  bind_rows() %>%
  arrange(order) %>%
  mutate(
    d_diff_2v1 =
      case_when(
        ## odds ratio to d
        test %in% "fact" & p < 0.05 ~ log(diff_2v1) * (sqrt(3) / pi),
        ## mean diff to d
        test %in% "cont" & p < 0.05 ~ compute.es::mes(m_2, m_1, sd_2, sd_1, .n_c123$sd_2, .n_c123$sd_1, verbose = FALSE)$d,  # nolint
        TRUE ~ NA_real_
      ),
    d_diff_3v1 =
      case_when(
        test %in% "fact" & p < 0.05 ~ log(diff_3v1) * (sqrt(3) / pi),
        test %in% "cont" & p < 0.05 ~ compute.es::mes(m_3, m_1, sd_3, sd_1, .n_c123$sd_3, .n_c123$sd_1, verbose = FALSE)$d,  # nolint
        TRUE ~ NA_real_
      ),
    d_diff_2v3 =
      case_when(
        test %in% "fact" & p < 0.05 ~ log(diff_2v3) * (sqrt(3) / pi),
        test %in% "cont" & p < 0.05 ~ compute.es::mes(m_2, m_3, sd_2, sd_3, .n_c123$sd_2, .n_c123$sd_3, verbose = FALSE)$d,  # nolint
        TRUE ~ NA_real_
      ),
    across(starts_with("d_diff_"), ~ round(.x, 2)),
    piff_2v1 = case_when(
      is.na(diff_2v1) ~ NA_character_,
      p < 0.05 & abs(x2_1) < bonf ~
        paste0("<b>", sprintf("%.2f", diff_2v1), "(", sprintf("%.2f", d_diff_2v1), ")", "</b>"),                             # nolint
      p < 0.05 & abs(x2_1) > bonf ~
        paste0(sprintf("%.2f", diff_2v1), "(", sprintf("%.2f", d_diff_2v1), ")"),                                            # nolint
      TRUE ~ NA_character_
    ),
    piff_3v1 = case_when(
      is.na(diff_3v1) ~ NA_character_,
      p < 0.05 & abs(x3_1) < bonf ~
        paste0("<b>", sprintf("%.2f", diff_3v1), "(", sprintf("%.2f", d_diff_3v1), ")", "</b>"),                             # nolint
      p < 0.05 & abs(x3_1) > bonf ~
        paste0(sprintf("%.2f", diff_3v1), "(", sprintf("%.2f", d_diff_3v1), ")"),                                            # nolint
      TRUE ~ NA_character_
    ),
    piff_2v3 = case_when(
      is.na(diff_2v3) ~ NA_character_,
      p < 0.05 & abs(x3_2) < bonf ~
        paste0("<b>", sprintf("%.2f", diff_2v3), "(", sprintf("%.2f", d_diff_2v3), ")", "</b>"),                             # nolint
      p < 0.05 & abs(x3_2) > bonf ~
        paste0(sprintf("%.2f", diff_2v3), "(", sprintf("%.2f", d_diff_2v3), ")"),                                            # nolint
      p > 0.05 & test %in% "fact" ~ NA_character_,
      TRUE ~ NA_character_
    ),
    across(starts_with("piff_"), ~ str_remove(.x, "NA(NA)")),
    across(
      starts_with("d_diff_"),
      ~ case_when(
        abs(.x) >= 0 & abs(.x) < 0.5 ~ "small",
        abs(.x) >= 0.5 & abs(.x) < 0.8 ~ "medium",
        abs(.x) >= 0.8 ~ "large"
      ),
      .names = "{.col}_fx"
    )
  )

### ============================================================================
### OUTCOMES
### ============================================================================
## DISTRIBUTION PLOTS OF ANOVA VARS
## -----------------------------------------------------------------------------
plot_iter <-
  lpa_t0t1_outc %>%
  dplyr::select(
    si, pb, tb, bhs, dep, anx, ac,
    str, age, sex, cov_sess
  ) %>%
  names() %>%
  purrr::set_names()

dens_fun <-
  function(dat, x, filly, facet) {
    ggplot(dat, aes(x = .data[[x]], fill = .data[[filly]])) + # nolint
      geom_density(alpha = 0.3) +                             # nolint
      facet_grid(.data[[facet]] ~ .) +                        # nolint
      theme_bw()                                              # nolint
  }

.plots_lpa_long <-
  map(plot_iter, ~ dens_fun(lpa_t0t1_outc, .x, "c", "final_session"))

plots_lpa_long <-
  wrap_plots(.plots_lpa_long, guides = "collect")

ggsave(
  "./figures/plots_lpa_long.png",
  plots_lpa_long
)

### ============================================================================
### CHANGE ANALYSIS STUFF
### ============================================================================
all_ch_plot <-
  lpa_t0t1_outc_w %>%
  select(c, si_ch, pb_ch, tb_ch, bhs_ch) %>%
  pivot_longer(cols = c(-c)) %>%
  mutate(
    name2 = toupper(name),
    name = case_when(
      name %in% "si_ch" ~ str_c("0", name),
      name %in% "pb_ch" ~ str_c("1", name),
      name %in% "tb_ch" ~ str_c("2", name),
      name %in% "bhs_ch" ~ str_c("3", name),
    )
  ) %>%
  ggplot(aes(x = value, fill = c)) +
  facet_grid(name2 ~ .) +
  geom_density() +
  theme_bw()

## DIFFERENCE CHANGE SCORE ANALYSIS
## =============================================================================
## https://www.theanalysisfactor.com/center-on-the-mean/
## https://www.theanalysisfactor.com/when-not-to-center-a-predictor-variable-in-regression/ # nolint
lpa_t0t1_outc_w <-
  lpa_t0t1_outc_w %>%
  mutate(across(
    c(pb_ch, tb_ch, bhs_ch),
    ~ scale(.x, scale = FALSE),
    .names = "{.col}_z"
  ))

lm_ch_class <-
  list(
    full = lpa_t0t1_outc_w,
    c1 = lpa_t0t1_outc_w %>% filter(c %in% 1),
    c2 = lpa_t0t1_outc_w %>% filter(c %in% 2),
    c3 = lpa_t0t1_outc_w %>% filter(c %in% 3)
    ) %>%
  map(~ .x %$% lm(si_ch ~ pb_ch * tb_ch * bhs_ch))

## DIFFERENCE CHANGE SCORE ANALYSIS WITH RANK TRANSFORMED SI_CH
## =============================================================================
lm_ch_r_class <- map(lm_ch_class, ~ update(.x, si_ch_r ~ .))

tt1 <- map(lm_ch_class, ~ update(.x, si_ch_r ~ pb_ch + tb_ch))
map(tt1, ~ summary(.x))

lm_ch_r_z_class <-
  map(lm_ch_class, ~ update(.x, si_ch_r ~ pb_ch_z * tb_ch_z * bhs_ch_z))

map(lm_ch_class, ~ summary(.x))
map(lm_ch_r_class, ~ summary(.x))
map(lm_ch_r_z_class, ~ summary(.x))

## diagnostic plots
plot(lm_ch_class$full)
plot(lm_ch_r_class$full)
plot(lm_ch_class$c1)
plot(lm_ch_r_class$c1)
plot(lm_ch_class$c2)
plot(lm_ch_r_class$c2)
plot(lm_ch_class$c3)
plot(lm_ch_r_class$c3)

## https://www.statology.org/how-to-calculate-the-p-value-of-an-f-statistic-in-r/ # nolint
## F statistic
f_mods_lm <-
  map(lm_ch_r_class, ~ as.data.frame(summary(.x)$fstatistic)) %>%
  map(
    ~ .x %>%
      rownames_to_column("F") %>%
      rename(val = 2)
  ) %>%
  bind_rows(.id = "c") %>%
  pivot_wider(
    names_from = F,
    values_from = val
  ) %>%
  mutate(
    p_value = pf(value, numdf, dendf, lower.tail = FALSE),
    p_value = case_when(
      p_value < 0.001 ~ "*p* < .001",
      p_value < 0.01 ~ "*p* < .01",
      p_value < 0.05 ~ "*p* < .05",
      TRUE ~ paste0("*p* = ", sprintf("%.3f", p_value))
    ),
    across(c(value), ~ sprintf("%.2f", .x)),
    F = paste0("*F*(", numdf, ", ", dendf, ") = ", value, ", ", p_value),
  )

r2_mods_lm <-
  list(
    r2 = map(lm_ch_r_class, ~ summary(.x)$r.squared),
    r2_adj = map(lm_ch_r_class, ~ summary(.x)$adj.r.squared)
  ) %>%
  bind_rows(.id = "type") %>%
  pivot_longer(-type) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(
    across(where(is.double), ~ str_replace(sprintf("%.2f", .x), "0\\.", "\\.")),
    r2_p = paste0("*R^2^* = ", r2),
    r2_adj_p = paste0("*R^2^~Adjusted~* = ", r2_adj),
    r2_print = paste0(r2_p, ", ", r2_adj_p)
  )

## PAIRED T-TEST
## -----------------------------------------------------------------------------
## TODO use MI data for this
.lpa_prepost_ttest <-
  list(
    full = list(
      si = pairwise.t.test.with.t.and.df(lpa_t0t1_outc$si, lpa_t0t1_outc$final_session, p.adjust.method = "bonf", paired = TRUE), # nolint
      pb = pairwise.t.test.with.t.and.df(lpa_t0t1_outc$pb, lpa_t0t1_outc$final_session, p.adjust.method = "bonf", paired = TRUE), # nolint
      tb = pairwise.t.test.with.t.and.df(lpa_t0t1_outc$tb, lpa_t0t1_outc$final_session, p.adjust.method = "bonf", paired = TRUE), # nolint
      bhs = pairwise.t.test.with.t.and.df(lpa_t0t1_outc$bhs, lpa_t0t1_outc$final_session, p.adjust.method = "bonf", paired = TRUE) # nolint
    ) %>%
      map(`[`, c("p.value", "t.value", "dfs")) %>%
      bind_rows(.id = "var") %>%
      mutate(across(everything(), ~ as.vector(.x))) %>%
      mutate(c = as.factor("0"), .after = var) %>%
      clean_names(),
    c123 = list(
      si = lpa_t0t1_outc %>%
        split(.$c) %>%
        map(~ pairwise.t.test.with.t.and.df(.x$si, .x$final_session, p.adjust.method = "bonf")), # nolint
      pb = lpa_t0t1_outc %>%
        split(.$c) %>%
        map(~ pairwise.t.test.with.t.and.df(.x$pb, .x$final_session, p.adjust.method = "bonf")), # nolint
      tb = lpa_t0t1_outc %>%
        split(.$c) %>%
        map(~ pairwise.t.test.with.t.and.df(.x$tb, .x$final_session, p.adjust.method = "bonf")), # nolint
      bhs = lpa_t0t1_outc %>%
        split(.$c) %>%
        map(~ pairwise.t.test.with.t.and.df(.x$bhs, .x$final_session, p.adjust.method = "bonf")) # nolint
    ) %>%
      unlist(recursive = FALSE) %>%
      map(`[`, c("p.value", "t.value", "dfs")) %>%
      bind_rows(.id = "var") %>%
      mutate(across(everything(), ~ as.vector(.x))) %>%
      separate(var, c("var", "c")) %>%
      mutate(c = as.factor(c)) %>%
      clean_names()
  ) %>%
  bind_rows(.id = ".drop")

.lpa_prepost_df <-
  list(
    full = lpa_t0t1_outc_w %>%
      select(c, ID, si_0, pb_0, tb_0, bhs_0, si_1, pb_1, tb_1, bhs_1) %>%
      ## group_by(c) %>%
      ## NOTE dep based means and standard deviations
      summarise(
        across(
          c(si_0, pb_0, tb_0, bhs_0, si_1, pb_1, tb_1, bhs_1),
          list(
            n = ~ sum(!is.na(.x)),
            m = ~ mean(.x, na.rm = TRUE),
            sd = ~ sd(.x, na.rm = TRUE)
          )
        ) %>%
          mutate(c = as.factor("0"), .before = 1)
      ),
    c123 = lpa_t0t1_outc_w %>%
      select(c, ID, si_0, pb_0, tb_0, bhs_0, si_1, pb_1, tb_1, bhs_1) %>%
      group_by(c) %>%
      ## NOTE class based means and standard deviations
      summarise(
        across(
          c(si_0, pb_0, tb_0, bhs_0, si_1, pb_1, tb_1, bhs_1),
          list(
            n = ~ sum(!is.na(.x)),
            m = ~ mean(.x, na.rm = TRUE),
            sd = ~ sd(.x, na.rm = TRUE)
          )
        )
      )
  ) %>%
  bind_rows(.id = ".drop")

.lpa_prepost_df <-
  list(
    full =
      lpa_t0t1_outc_w %>%
        select(c, ID, si_0, pb_0, tb_0, bhs_0, si_1, pb_1, tb_1, bhs_1) %>%
        pivot_longer(-c(c, ID), names_to = c("var", "time"), names_sep = "_") %>% # nolint
        pivot_wider(names_from = time, values_from = value, names_glue = "t{.name}") %>% # nolint
        split(.$var) %>%
        map(~ .x %$% DescTools::MeanDiffCI(t0, t1, na.rm = TRUE)) %>%
        bind_rows(.id = "var") %>%
        mutate(c = "full"),
    c123 =
      lpa_t0t1_outc_w %>%
        select(c, ID, si_0, pb_0, tb_0, bhs_0, si_1, pb_1, tb_1, bhs_1) %>%
        pivot_longer(-c(c, ID), names_to = c("var", "time"), names_sep = "_") %>% # nolint
        pivot_wider(names_from = time, values_from = value, names_glue = "t{.name}") %>% # nolint
        split(.$var) %>%
        map(
          ~ .x %>%
            split(.$c) %>%
            map(~ .x %$% DescTools::MeanDiffCI(t0, t1, na.rm = TRUE)) %>%
            bind_rows(.id = "c")
        ) %>%
        bind_rows(.id = "var")
  ) %>%
  bind_rows(.id = ".drop")

.lpa_prepost_desc <-
  .lpa_prepost_df %>%
  pivot_longer(-c(c, .drop)) %>%
  mutate(
    time = str_extract(name, "\\d"),
    type = str_extract(name, "(n|sd|m)$"),
    var = str_extract(name, "^(si|pb|tb|bhs)"),
    across(
      value:time,
      ~ as.numeric(.x)
    )
  ) %>%
  select(-name) %>%
  pivot_wider(
    names_from = c(type),
    values_from = value
  ) %>%
  mutate(order = case_when(
    var %in% "si" ~ 1,
    var %in% "pb" ~ 2,
    var %in% "tb" ~ 3,
    var %in% "bhs" ~ 4
  ), .before = 1) %>%
  arrange(order, time)

.lpa_prepost_fx_mes <-
  pmap(
    list(
      m_1 = paste0(
        ".lpa_prepost_df$",
        c("si", "pb", "tb", "bhs"),
        "_",
        rep(0, 4), "_m"
      ) %>%
        map(~ eval(parse(text = .x))),
      m_2 = paste0(
        ".lpa_prepost_df$",
        c("si", "pb", "tb", "bhs"),
        "_",
        rep(1, 4), "_m"
      ) %>%
        map(~ eval(parse(text = .x))),
      sd_1 = paste0(
        ".lpa_prepost_df$",
        c("si", "pb", "tb", "bhs"),
        "_",
        rep(0, 4), "_sd"
      ) %>%
        map(~ eval(parse(text = .x))),
      sd_2 = paste0(
        ".lpa_prepost_df$",
        c("si", "pb", "tb", "bhs"),
        "_",
        rep(1, 4), "_sd"
      ) %>%
        map(~ eval(parse(text = .x))),
      n_1 = paste0(
        ".lpa_prepost_df$",
        c("si", "pb", "tb", "bhs"),
        "_",
        rep(0, 4), "_n"
      ) %>%
        map(~ eval(parse(text = .x))),
      n_2 = paste0(
        ".lpa_prepost_df$",
        c("si", "pb", "tb", "bhs"),
        "_",
        rep(1, 4), "_n"
      ) %>%
        map(~ eval(parse(text = .x)))
    ),
    .f = function(m_1, m_2, sd_1, sd_2, n_1, n_2) {
      compute.es::mes(m_1, m_2, sd_1, sd_2, n_1, n_2)
    }
  )

.lpa_prepost_d <-
  .lpa_prepost_fx_mes %>%
  set_names(c("si", "pb", "tb", "bhs")) %>%
  map(
    ~ .x %>%
      rowid_to_column("c") %>%
      mutate(c = as.factor(c - 1)) %>%
      select(c, d, l.d, u.d)
  ) %>%
  bind_rows(.id = "var") %>%
  mutate(
    across(
      d,
      ~ case_when(
        abs(.x) >= 0 & abs(.x) < 0.5 ~ "small",
        abs(.x) >= 0.5 & abs(.x) < 0.8 ~ "medium",
        abs(.x) >= 0.8 ~ "large"
      ),
      .names = "{.col}_fx"
    )
  )

.lpa_prepost_kbl <-
  .lpa_prepost_desc %>%
  full_join(.lpa_prepost_ttest) %>%
  full_join(.lpa_prepost_d) %>%
  select(-.drop)

lpa_prepost_kbl <-
  .lpa_prepost_kbl %>%
  group_by(var, c) %>%
  mutate(
    .m_d = sprintf("%.2f", m[2] - m[1]), .before = d
  ) %>%
  ungroup() %>%
  mutate(
    across(
      c(m, sd, t_value, d),
      ~ sprintf("%.2f", .x)
    ),
    m_sd = paste0(m, "(", sd, ")"),
    t_df = case_when(
      p_value < 0.001 ~ paste0(t_value, "(", dfs, ")", "<sup>***</sup>"), # nolint
      p_value < 0.01 ~  paste0(t_value, "(", dfs, ")", "<sup>**</sup>"), # nolint
      p_value < 0.05 ~  paste0(t_value, "(", dfs, ")", "<sup>*</sup>"), # nolint
      p_value > 0.05 ~  paste0(t_value, "(", dfs, ")", ""),
      TRUE ~ NA_character_
    ),
    m_d = paste0(.m_d, "(", d, ")"),
    var_t = case_when(
      time %in% 0 ~ paste0(toupper(var), "<sub>1</sub>"),
      time %in% 1 ~ paste0(toupper(var), "<sub>2</sub>")
    ), .after = var
  ) %>%
  select(order, c, time, var, var_t, n, m_sd, t_df, m_d, d_fx) %>% # nolint
  pivot_wider(
    names_from = c,
    names_glue = "c{c}_{.value}",
    values_from = n:last_col()
  ) %>%
  select(order, time, var, var_t, starts_with(c("c0", "c1", "c2", "c3")))

lpa_pp_long <-
  .lpa_prepost_kbl %>%
  mutate(
    across(
      c(m, sd, t_value, d),
      ~ sprintf("%.2f", .x)
    ),
    m_sd = paste0(m, "(", sd, ")"),
    p_print = case_when(
      round(p_value, 2) <= 0.05 ~ "*p* < .05",
      round(p_value, 2) <= 0.01 ~ "*p* < .01",
      round(p_value, 3) <= 0.001 ~ "*p* < .001",
      TRUE ~ paste0("*p* = ", str_remove(sprintf("%.2f", p_value), "^0"))
    ),
    t_df = paste0("*t*(", dfs, ")", " = ", t_value, ", ", p_print),
  ) %>%
  select(c, time, m, sd, var, t_df, d_fx)

### ============================================================================
### ABSTRACT DF
### ============================================================================
## NOTE this is okay
abst_df <-
  lpa_t0_og %>%
  summarise(
    N = sum(!is.na(c)),
    age_m = round(mean(age, na.rm = TRUE), 2),
    age_sd = round(sd(age, na.rm = TRUE), 2),
    prc_1 = round(sum(!is.na(c[c %in% 1])) / N * 100, 2),
    prc_2 = round(sum(!is.na(c[c %in% 2])) / N * 100, 2),
    prc_3 = round(sum(!is.na(c[c %in% 3])) / N * 100, 2),
  )

## CORRELATION TABLES
## =============================================================================
corr_des <-
  lpa_t0t1_outc_w %>%
  select(
    !ends_with(c("_z", "_r", "_ch")) & where(is.numeric),
    -session_0, sex, starts_with("hx"), pastpsyc, atsi, cald
  ) %>%
  rename_with(~ str_remove(.x, "^isas_")) %>%
  rename(
    session_0 = session_1,
    DSHage_0 = `3a`,
    DSHn_0 = dsh_n,
    pastpsyc_0 = pastpsyc,
    atsi_0 = atsi,
    cald_0 = cald
  ) %>%
  rename_with(~ str_c(.x, "_0"), starts_with(c("hx", "age", "sex", "ac", "soc", "int"))) %>% # nolint
  mutate(across(everything(), as.numeric))

corr_tbl <-
  corstars(corr_des, result = "html") %>%
  rownames_to_column("var") %>%
  separate(var, c("var", "time")) %>%
  mutate(
    time = as.integer(time),
    var = case_when(
      var %in% c("pb", "ac", "tb", "dep", "anx", "str", "bhs") ~ paste0(toupper(var), "<sub>", time + 1, "</sub>"), # nolint
      var %in% "si" ~ paste0("MSSI", "<sub>", time + 1, "</sub>"),
      var %in% "sex" ~ "Sex<sup>a</sup>",
      var %in% "pastpsyc" ~ "Treatment<sub>Hx</sub><sup>b</sup>",
      var %in% "hxtrauma" ~ "Trauma<sub>Hx</sub><sup>b</sup>",
      var %in% "hxsuicide" ~ "Suicide<sub>Hx</sub><sup>b</sup>",
      var %in% "hxdsh" ~ "DSH<sub>Hx</sub><sup>b</sup>",
      var %in% "hxsa" ~ "SA<sub>Hx</sub><sup>b</sup>",
      var %in% "atsi" ~ "ATSI<sup>b</sup>",
      var %in% "cald" ~ "CALD<sup>b</sup>",
      var %in% "DSHn" ~ "DSH<sub>n</sub>",
      var %in% "session" ~ "Session<sub>n</sub>",
      var %in% "social" ~ "ISAS<sub>inter</sub>",
      var %in% "intra" ~ "ISAS<sub>intra</sub>",
      var %in% "age" ~ "Age",
      var %in% "DSHage" ~ "DSH<sub>age</sub>",
      TRUE ~ var
    )
  ) %>%
  select(-time) %>%
  rename(Variable = var)

names(corr_tbl) <-
  corr_tbl %>%
  select(Variable) %>%
  t() %>%
  as.vector() %>%
  append(., "Variable", 0)

## =============================================================================
corr_des %>%
  cor.ci() %>%
  cor.plot.upperLowerCi()

### ============================================================================
### PLOTS AND GRAPHS
### ============================================================================
.lpa_plot <-
  mplus_list["c3_2_lpa_si2_v1.out"]

lpa_fig_plot <-
  plotMixtures(
    .lpa_plot,
    coefficients = "stdyx.standardized",
    ) +
  theme_ipsum(
    axis_title_face = "bold",
    axis_title_size = 12,
    axis_title_just = "c"
  ) +
  labs(
    x = "ITS Factors",
    y = "Standardised Paramater Estimates"
  ) +
  scale_x_discrete(label = c("AC", "BHS", "PB", "TB"))

ggsave(
  "./figures/lpa_fig_plot.png",
  lpa_fig_plot
)

### ============================================================================
### PARTCIPANTS FOR METHODS SECTION
### ============================================================================
## TODO update this for current data set
## consenters
.phd_sp_full <-
  read_csv(here(.phd_sp_full_df)) %>%
  filter(sp_episode %in% 1)

.lpa_consenters <-
  .phd_sp_full %>%
  group_by(uci) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(
    final_session = if_else(date %in% max(date) & days_since_session > 90, 1, 0), # nolint
    .after = session
  ) %>%
  ungroup() %>%
  mutate(
    ## calculate percent missing on mssi screener items
    mssi_scrn_miss2 = rowMeans(is.na(across(mssi_1:mssi_4))),
    ## calculate whether below cut-off (co) no %in% 0, yes = 1
    mssi_co = case_when(
      mssi_1 %in% (c(1, 0, NA)) &
        mssi_2 %in% (c(1, 0, NA)) &
          mssi_3 %in% (c(0, NA))    &
          mssi_4 %in% (c(0, NA))    ~ 1,
      TRUE ~ 0
    ),
    ## if items 5--18 == NA and mssi_co below cut-off and have no screening
    ## items missing, replace with 0, else treat as NA
    across(mssi_5:mssi_18, ~ case_when(
      is.na(.) & mssi_scrn_miss2 < 0.25 & mssi_co %in% 1 ~ 0,
      TRUE ~ as.numeric(.)
    ), .names = "rc_{.col}")
  ) %>%
  .pro_mean("si", c(mssi_1:mssi_4, starts_with("rc_mssi"))) %>%
  select(ID, session, final_session, sp_episode, consent, si)

lpa_consenters <-
  .lpa_consenters %>%
  filter(si > 0, final_session %in% 0, session %in% 1) %>%
  group_by(consent) %>%
  summarise(n = n()) %>%
  filter(!consent %in% 1)

## SEX AND AGE DIFFERENCE
## -----------------------------------------------------------------------------
## TODO keep
## sex diff with age
.lpa_sex_diff <-
  lpa_fin %>%
  filter(!is.na(c)) %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  select(ID, c, sex, age)
## TODO keep
sltb_age_sex_t <-
  pairwise.t.test.with.t.and.df(
    .lpa_sex_diff$age,
    .lpa_sex_diff$sex,
    p.adjust.method = "bonf"
  )

lpa_sex_age_diff <-
  .lpa_sex_diff %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    across(
      age,
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(.x, na.rm = TRUE)
      ))
  ) %>%
  filter(sex %in% c(0:2)) %>%
  mutate(
    across(c(age_mean, age_median, age_sd), ~ sprintf("%.2f", .x))
  )

### ============================================================================
### WRITE IMAGE AND DATA
### ============================================================================
save.image("./data/proc/lpa_rda.RData")
saveRDS(mget(ls()), "./data/proc/lpa_df.rds")

### ============================================================================
### END SCRIPT
### ============================================================================

## CATEGORICAL CONTRASTS WITH CIs
## -----------------------------------------------------------------------------

.fu_multnom <-
  .fu %>%
  select(sex, ends_with("2"), c) %>%
  mutate(
    across(where(is.factor), ~ tolower(as.character(.x))),
    sex = if_else(sex %in% "female", "yfemale", sex),
  ) %>%
  pivot_longer(-c)

.fu_or_multinom <-
  list(
    x1 =
      .fu_multnom %>%
        split(.$name) %>%
        map(~ .x %>% mutate(c = relevel(as.factor(c), ref = "1"))),
    x2 =
      .fu_multnom %>%
        split(.$name) %>%
        map(~ .x %>% mutate(c = relevel(as.factor(c), ref = "2"))),
    x3 =
      .fu_multnom %>%
        split(.$name) %>%
        map(~ .x %>% mutate(c = relevel(as.factor(c), ref = "3")))
  ) %>%
  map(
    ~ .x %>%
      map(
        ~ .x %$%
          nnet::multinom(c ~ value) %>%
          tidy(exponentiate = TRUE, conf.int = TRUE) %>%
          filter(!term %in% "(Intercept)") %>%
          select(y.level:estimate, p.value:conf.high) %>%
          mutate(
            term = str_remove(term, "^value"),
            p.value = 3 * p.value
          ) %>%
          rename(c = 1, or = estimate)
      ) %>%
      bind_rows(.id = "var")
  ) %>%
  bind_rows(.id = "cont") %>%
  filter(!cont %in% "x1") %>%
  unite(cont, c(cont, c)) %>%
  filter(!cont %in% "x3_2") %>%
  rename(
    mdif_or = or,
    ci_l = conf.low,
    ci_u = conf.high,
    p_or = p.value
  ) %>%
  pivot_wider(
    names_from = cont,
    values_from = mdif_or:last_col()
  )

## NOTE the below can be calculated form the fu_test table. Do this, and remove
## the latter as superfulous
.fu_chisq_bonf_ds <-
  list(
    x2_1 = .chiq_bonf_21 %>%
      mutate(across(where(is.character), ~ as.factor(.x))) %>%
      pivot_longer(-c) %>%
      split(.$name) %>%
      map(~ .x %$% tidy(chisq.test(value, c))) %>%
      unlist(recursive = FALSE) %>%
      as_tibble() %>%
      select(ends_with(c("statistic"))) %>%
      pivot_longer(everything(), names_to = c("var", "stat"), names_sep = "\\.") %>% # nolint
      pivot_wider(names_from = stat, values_from = value) %>%
      split(.$var) %>%
      map(
        ~ .x %$%
          compute.es::chies(statistic, nrow(.chiq_bonf_21)) %>%
          select(d, l.d, u.d)
      ) %>%
      bind_rows(.id = "var"),
    x3_1 = .chiq_bonf_31 %>%
      mutate(across(where(is.character), ~ as.factor(.x))) %>%
      pivot_longer(-c) %>%
      split(.$name) %>%
      map(~ .x %$% tidy(chisq.test(value, c))) %>%
      unlist(recursive = FALSE) %>%
      as_tibble() %>%
      select(ends_with(c("statistic"))) %>%
      pivot_longer(everything(), names_to = c("var", "stat"), names_sep = "\\.") %>% # nolint
      pivot_wider(names_from = stat, values_from = value) %>%
      split(.$var) %>%
      map(
        ~ .x %$%
          compute.es::chies(statistic, nrow(.chiq_bonf_31)) %>%
          select(d, l.d, u.d)
      ) %>%
      bind_rows(.id = "var"),
    x2_3 = .chiq_bonf_23 %>%
      mutate(across(where(is.character), ~ as.factor(.x))) %>%
      pivot_longer(-c) %>%
      split(.$name) %>%
      map(~ .x %$% tidy(chisq.test(value, c))) %>%
      unlist(recursive = FALSE) %>%
      as_tibble() %>%
      select(ends_with(c("statistic"))) %>%
      pivot_longer(everything(), names_to = c("var", "stat"), names_sep = "\\.") %>% # nolint
      pivot_wider(names_from = stat, values_from = value) %>%
      split(.$var) %>%
      map(
        ~ .x %$%
          compute.es::chies(statistic, nrow(.chiq_bonf_23)) %>%
          select(d, l.d, u.d)
      ) %>%
      bind_rows(.id = "var")
  ) %>%
  bind_rows(.id = "cont") %>%
  filter(!var %in% "planned_final2") %>%
  rename(d_l = l.d, d_u = u.d) %>%
  pivot_wider(names_from = cont, values_from = d:d_u) %>%
  select(var, ends_with(c("2_1", "3_1", "2_3")))

.fu_chisq_bonf_ds %>%
  select(starts_with(c("d_u_", "d_l_"))) %>%
  rename_with(
    ~ str_replace_all(.x, c("(d)(_u|_l)(_x\\d_\\d)" = "\1\3\2")),
    everything()
    ## starts_with(c("d_u_", "d_l_"))
)

.fu_d_fact <-
  full_join(.fu_or_multinom, .fu_chisq_bonf_ds) %>%
  mutate(test = "fact")

## CONTINUOUS CONTRASTS WITH CIs
## -----------------------------------------------------------------------------

.lpa_msd_full <-
  .lpa_anova %>%
  map(
    ~ .x %>%
      split(.$name) %>%
      map(
        ~ .x %>%
          group_by(c) %>%
          summarise(
            n = n(),
            across(
              value,
              list(
                m = ~ mean(.x, na.rm = TRUE),
                sd = ~ sd(.x, na.rm = TRUE)
              )
            )
          ) %>%
          rename_with(~ str_remove(.x, "value_")) %>%
          pivot_wider(
            names_from = c,
            values_from = n:last_col()
          )
      ) %>%
      bind_rows(.id = "var") %>%
      mutate(
        d_x2_1 = compute.es::mes(m_2, m_1, sd_2, sd_1, n_2, n_1, verbose = FALSE)$d, # nolint
        d_x3_1 = compute.es::mes(m_3, m_1, sd_3, sd_1, n_3, n_1, verbose = FALSE)$d, # nolint
        d_x2_3 = compute.es::mes(m_2, m_3, sd_2, sd_3, n_2, n_3, verbose = FALSE)$d, # nolint
        d_x2_1_l = compute.es::mes(m_2, m_1, sd_2, sd_1, n_2, n_1, verbose = FALSE)$l.d, # nolint
        d_x3_1_l = compute.es::mes(m_3, m_1, sd_3, sd_1, n_3, n_1, verbose = FALSE)$l.d, # nolint
        d_x2_3_l = compute.es::mes(m_2, m_3, sd_2, sd_3, n_2, n_3, verbose = FALSE)$l.d, # nolint
        d_x2_1_u = compute.es::mes(m_2, m_1, sd_2, sd_1, n_2, n_1, verbose = FALSE)$u.d, # nolint
        d_x3_1_u = compute.es::mes(m_3, m_1, sd_3, sd_1, n_3, n_1, verbose = FALSE)$u.d, # nolint
        d_x2_3_u = compute.es::mes(m_2, m_3, sd_2, sd_3, n_2, n_3, verbose = FALSE)$u.d, # nolint
      )
  ) %>%
  set_names(c("t0", "t1", "ch")) %>%
  bind_rows(.id = "time") %>%
  filter(time %in% "t0") %>%
  select(-time)

.mean_diff_ci <-
  list(
    x2_1 =
      .lpa_anova %>%
        map(
          ~ .x %>%
            mutate(c = as.character(c)) %>%
            filter(!c %in% 3) %>%
            pivot_wider(
              names_from = c, values_from = value,
              names_glue = "c{.name}"
            ) %>%
            split(.$name) %>%
            map(
              ~ .x %$%
                DescTools::MeanDiffCI(c2, c1, na.rm = TRUE)
            ) %>%
            bind_rows(.id = "var")
        ) %>%
        set_names(c("t0", "t1", "ch")),
    x3_1 =
      .lpa_anova %>%
        map(
          ~ .x %>%
            mutate(c = as.character(c)) %>%
            filter(!c %in% 2) %>%
            arrange(desc(c)) %>%
            pivot_wider(
              names_from = c, values_from = value,
              names_glue = "c{.name}"
            ) %>%
            split(.$name) %>%
            map(
              ~ .x %$%
                DescTools::MeanDiffCI(c3, c1, na.rm = TRUE)
            ) %>%
            bind_rows(.id = "var")
        ) %>%
        set_names(c("t0", "t1", "ch")),
    x2_3 =
      .lpa_anova %>%
        map(
          ~ .x %>%
            mutate(c = as.character(c)) %>%
            filter(!c %in% 1) %>%
            arrange(c) %>%
            pivot_wider(
              names_from = c, values_from = value,
              names_glue = "c{.name}"
            ) %>%
            split(.$name) %>%
            map(
              ~ .x %$%
                DescTools::MeanDiffCI(c2, c3, na.rm = TRUE)
            ) %>%
            bind_rows(.id = "var")
        ) %>%
        set_names(c("t0", "t1", "ch"))
  ) %>%
  map(~ .x %>% bind_rows(.id = "time")) %>%
  bind_rows(.id = "cont") %>%
  filter(time %in% "t0") %>%
  select(-time) %>%
  rename(
    mdif_or = meandiff,
    ci_l = lwr.ci,
    ci_u = upr.ci
  ) %>%
  pivot_wider(
    names_from = cont,
    values_from = mdif_or:ci_u,
  )

.fu_d_cont <-
  full_join(.lpa_msd_full, .mean_diff_ci) %>%
  mutate(test = "cont")

## MERGE ABOVE
## -----------------------------------------------------------------------------

.fu_d_df <-
  full_join(.fu_d_cont, .fu_d_fact) %>%
  rename(var3 = var, value3 = term) %>%
  select(-c(n_1:sd_3)) %>%
  filter(!var3 %in% "planned_final2") %>%
  mutate(
    var3 = str_remove(var3, "2$"),
    var3 = if_else(var3 %in% "session", "total_sessions", var3)
  )


.fff <-
  full_join(.fu_test_full, .fu_d_df)

.fu_d_df_2 <-
  .fu_d_df %>%
  mutate(
    across(
      where(is.numeric) & !starts_with("p_"),
      ~ sprintf("%.2f", .x)
    )
  )
