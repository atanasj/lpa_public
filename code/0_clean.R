### ============================================================================
### SOURCE SCRIPTS
### ============================================================================
rm(list = ls(all = TRUE))
source("./code/func.R")
source("./code/load.R")

### ============================================================================
### NOTES ON DATA CLEANSING
### ============================================================================

### ============================================================================
### LOAD DATA
### ============================================================================

.phd_sp_full <-
  read_csv(here(.phd_sp_full_df))

## create list from uci that have a first session recorded
.phd_target <-
  ## get data
  .phd_sp_full %>%
  ## filter
  filter(
    ## first episode only
    sp_episode == 1,
    ## must provide consent
    consent == 1,
    ## must have admission measures
    session == 1
  ) %>%
  select(uci) %>%
  ## make the list
  deframe()

## create df with at least first session
.phd_sp <-
  .phd_sp_full %>%
  ## filter(sp_episode == 1,
  ##        consent == 1,
  ##        ## uci %in% .phd_target
  ##        ) %>%
  arrange(uci, date, session)

## .phd_sp_full %>%
.phd_sp %>%
  filter(if_all(c(sp_episode, session), ~.x %in% 1)) %>%
  group_by(consent) %>%
  summarise(n = n())

## inspect lpa cases by number of sessions
.phd_sp_full %>%
## .phd_sp %>%
  filter(session == 1 | final_session == 1) %>%
  group_by(uci) %>%
  add_count(name = "n_sess") %>%
  ungroup() %>%
  group_by(n_sess) %>%
  summarise(n_obs = n_distinct(uci))

### ============================================================================
### TRANSFORM VARIABLES IN PREPARATION FOR ANALYSES
### ============================================================================

## participants with at least two time points / sessions, and potentially only
## participants with "clinical caseness" e.g., [@saunders2016]
## NOTE see GPH analysis for long to wide
## NOTE see cleansing project for other recoding of INQ
.phd_sp <-
  .phd_sp %>%
  mutate(
    ## calculate percent missing on mssi screener items
    mssi_scrn_miss2 = rowMeans(is.na(across(mssi_1:mssi_4))),
    ## calculate whether below cut-off (co) no == 0, yes = 1
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
      is.na(.) & mssi_scrn_miss2 < 0.25 & mssi_co == 1 ~ 0,
      TRUE ~ as.numeric(.)
    ), .names = "rc_{.col}")
  ) %>%
  ## recode and create new bhs variables
  ## NOTE this was intended for use with SEM, but was not used in final anslysis
  mutate(
    across(c(bhsSF_1r, bhsSF_2, bhsSF_3, bhsSF_4r), ~ case_when(
      . <= 3 ~ 0,
      . >= 4 ~ 1,
      TRUE ~ NA_real_
    ),
    .names = "rc_{.col}"
    ),
    bhsrs_6 = coalesce(bhs_6r, rc_bhsSF_1r),
    bhsrs_7 = coalesce(bhs_7r, rc_bhsSF_2),
    bhsrs_9 = coalesce(bhs_9r, rc_bhsSF_3),
    bhsrs_15 = coalesce(bhs_15r, rc_bhsSF_4r),
    ## BHS-20 positively worded items
    across(
      num_range("bhs_", c(1, 3, 5, 6, 8, 10, 13, 15, 19)),
      ~ case_when(
        .x %in% 1 ~ 0,
        .x %in% 2 ~ 1,
        TRUE ~ NA_real_
      ),
      .names = "{.col}r"
    ),
    ## BHS-20 negatively worded items
    across(
      num_range("bhs_", c(2, 4, 7, 9, 11, 12, 14, 16, 17, 18, 20)),
      ~ case_when(
        .x %in% 2 ~ 0,
        .x %in% 1 ~ 1,
        TRUE ~ NA_real_
      ),
      .names = "{.col}r"
    )
  ) %>%
  group_by(uci) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup() %>%
  arrange(ID, uci, session) %>%
  select(uci, ID, everything())

.phd_sp %>%
  group_by(uci) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup() %>%
  arrange(ID, uci, session) %>%
  filter(session < 6) %>%
  select(uci, ID, session, contains("mssi")) %>%
  .pro_mean("si", c(mssi_1:mssi_4, starts_with("rc_mssi"))) %>%
  select(session, si) %>%
  split(.$session) %>%
  map(
    ~ .x %>%
      summary()
  )

phd_sp <-
  ## NOTE this was used when calculating the means the LPA
  ## phd_sp_means %>%
  .phd_sp %>%
  group_by(uci) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup() %>%
  arrange(ID, uci, session) %>%
  select(uci, ID, everything()) %>%
  group_by(ID) %>%
  mutate(
    ## across(c(c:sex, atsi:closehxdsh, suicide_attempts), as.factor),
    ## create final session variable
    final_session = if_else(date == max(date) & days_since_session > 90, 1, 0),
    .after = session
  ) %>%
  ungroup() %>%
  mutate(across(
    c(bhsSF_1r, bhsSF_2, bhsSF_3, bhsSF_4r),
    .fns = function(x) {
      x - 1
    },
    .names = "{.col}_lessone"
  )) %>%
  ## get bhs sum of protocol 1
  .pro_mean("bhs_m1", starts_with("bhs_") & ends_with("r")) %>%
  ## get bhs sum of protocol 2
  .pro_mean("bhs_m2", starts_with("bhsSF_") & ends_with("lessone")) %>%
  .pro_mean("pb", inq_1:inq_6) %>%
  .pro_mean("tb", c(inq_7r:inq_15r, inq_9, inq_11:inq_12)) %>%
  ## caolesce protocol 1 and 2 into one score
  mutate(bhs = coalesce(bhs_m1, bhs_m2)) %>%
  .pro_mean("ac", c(acss_8r:acss_13r, acss_7, acss_11, acss_14, acss_19)) %>%
  ## create means of appropriate vars
  .pro_mean("dep", num_range("dass_", c(3, 5, 10, 13, 16, 17, 21))) %>%
  .pro_mean("anx", num_range("dass_", c(2, 4, 7, 9, 15, 19, 20))) %>%
  .pro_mean("str", num_range("dass_", c(1, 6, 8, 11, 12, 14, 18))) %>%
  .pro_mean("si", c(mssi_1:mssi_4, starts_with("rc_mssi"))) %>%
  .pro_mean("is_aff_reg", num_range("isas_8.", c(1, 14, 27))) %>%
  ## isas interpersonal bonding
  .pro_mean("is_int_bnd", num_range("isas_8.", c(2, 15, 28))) %>%
  ## isas self-punishment
  .pro_mean("is_slf_pun", num_range("isas_8.", c(3, 16, 29))) %>%
  ## isas self-care
  .pro_mean("is_slf_car", num_range("isas_8.", c(4, 17, 30))) %>%
  ## isas anti-dissociation
  .pro_mean("is_ant_dis", num_range("isas_8.", c(5, 18, 31))) %>%
  ## isas anti-suicide
  .pro_mean("is_ant_sui", num_range("isas_8.", c(6, 19, 32))) %>%
  ## isas sensation-seeking
  .pro_mean("is_sen_see", num_range("isas_8.", c(7, 20, 33))) %>%
  ## isas peer-bonding
  .pro_mean("is_peer_bnd", num_range("isas_8.", c(8, 21, 34))) %>%
  ## isas interpersonal influence
  .pro_mean("is_int_inf", num_range("isas_8.", c(9, 22, 35))) %>%
  ## isas toughness
  .pro_mean("is_tough", num_range("isas_8.", c(10, 23, 36))) %>%
  ## isas marking distress
  .pro_mean("is_mrk_dis", num_range("isas_8.", c(11, 24, 37))) %>%
  ## isas revenge
  .pro_mean("is_revenge", num_range("isas_8.", c(12, 25, 38))) %>%
  ## isas autonomy
  .pro_mean("is_autnmy", num_range("isas_8.", c(13, 36, 39))) %>%
  ## isas intrapersonal-factor
  .pro_mean("isas_intra", num_range(
    "isas_8.",
    ## intrapersonal-influence
    c(
      1, 3, 5, 6, 11, 14, 16, 18,
      19, 24, 27, 29, 31, 32, 37
    )
  )) %>%
  ## isas social-factors
  .pro_mean("isas_social", num_range(
    "isas_8.",
    ## social-interpersonal
    c(
      2, 4, 7, 8, 9, 10, 12, 13, 15, 17,
      20, 21, 22, 23, 25, 26, 28, 30,
      33, 34, 35, 36, 38, 39
    )
  ))

phd_sp %>%
  saveRDS(here(.phd_sp_df))

phd_sp %>%
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


lpa_fit <-
  phd_sp %>%
  filter(session == 1) %>%
  select(
    ID, sex, age, atsi:suicide_attempts, protocol,
    pb:isas_social, everything()
  )

lpa_fit %>%
  saveRDS(here(.lpa_fit_df))

lta_fit <-
  phd_sp %>%
  filter(session == 1 | final_session == 1) %>%
  mutate(
    drop = if_else(
      session == 1 & final_session == 1,
      1, 0
    ), .after = ID
  ) %>%
  ## select those with more than one session
  filter(drop == 0) %>%
  group_by(ID) %>%
  ## drop lone ranger ID 285
  filter(max(session) > 1) %>%
  ungroup() %>%
  select(
    uci, ID, sex, age, atsi:suicide_attempts, protocol,
    pb:isas_social, everything()
  )

lta_fit %>%
  saveRDS(here(.lta_fit_df))

### ============================================================================
### END SCRIPT: clean.R
### ============================================================================

show <-
  lpa_fit %>%
  filter(si > 0) %>%
  filter(final_session == 0) %>%
  select(si) %>%
  glimpse()
  group_by(consent) %>%
  summarise(n=n())
