## MULTIPLE IMPUTATION
## =============================================================================
## need to have the same data set
rm(list = ls(all = TRUE))
source("./code/func.R")
source("./code/load.R")
## load("./data/proc/lpa_rda.RData")
load.image("./data/proc/mi_rda.RData")
library(mice)
library(miceadds)           # for anova analysis on mids obejcts
library(broom.mixed)

lpa_t0t1 <-
  lpa_fin %>%
  ## filter no class var
  filter(!is.na(c)) %>%
  group_by(ID) %>%
  ## filter first and last session
  filter(session %in% 1 | final_session %in% 1) %>%
  ungroup() %>%
  mutate(
    ## replace NA session value
    session = if_else(is.na(session), -999, session),
    ## create conditional variable to drop those with only one session
    drop = if_else(session %in% 1 & final_session %in% 1, 1, 0), .after = ID,
    across(c(isas_1i, isas_1m), ~ str_extract(.x, "\\d*")),
    across(starts_with("isas_1"), as.numeric),
    ## prorate isas higher-order factors for comparison
    isas_intra = isas_intra / 5,
    isas_social = isas_social / 8
  ) %>%
  rowwise() %>%
  mutate(
    ## count the number of methods of dsh
    dsh_n = sum(c_across(starts_with("isas_1")) > 0, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  ## select those with more than one session
  filter(drop %in% 0) %>%
  ## select vars for mi
  select(
    c, ID, age, sex, session, final_session,
    planned_final, pb, ac, bhs, tb, dep, anx, si,
    isas_intra, isas_social, atsi, cald, hxtrauma, hxsuicide, hxdsh,
    isas_3a, # age onset
    dsh_n
  ) %>%
  group_by(ID) %>%
  mutate(
    age = min(age),
    planned_final = max(planned_final)
  ) %>%
  ungroup()

.lpa_t0t1_mi <-
  lpa_t0t1 %>%
  panelr::panel_data(
            id = ID,
            wave = final_session,
            ) %>%
  panelr::widen_panel(
            separator = "_",
            ignore.attributes = FALSE,
            varying = NULL) %>%
  ungroup() %>%
  ## need to convert to data frame as panelr creates a list object
  as.data.frame()

lpa_t0t1_mi <-
  .lpa_t0t1_mi %>%
  select(-ac_1) %>%
  rename(ac = ac_0) %>%
  ## use age at intake, drop secon ac as not routinely collected
  mutate(session_1 = if_else(session_1 %in% -999, NA_real_, session_1))

lpa_t0t1_mi %>% str()
lpa_t0t1_mi %>% summary()

lpa_t0t1_mi <-
  lpa_t0t1_mi %>%
  mutate(
    across(
      c(
        c, ID, sex, session_0, session_1,
        planned_final, atsi, cald, hxtrauma, hxsuicide, hxdsh,
        isas_3a
      ),
      ~ as.integer(as.character(.x))
    )
  )

imp_start <- mice(lpa_t0t1_mi, maxit = 0)

meth <- imp_start$method
pred <- imp_start$predictorMatrix
## pred[, colnames(pred) %in% c("ID")] <- 0
pred[, colnames(pred) %in% c("ID", "c")] <- 0
pred


imp_t0t1_5 <-
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
  imp_t0t1_5 %>%
  complete(., action = "long", include = TRUE) %>%
  group_by(ID) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  mutate(
    si_ch = si_0 - si_1,
    pb_ch = pb_0 - pb_1,
    tb_ch = tb_0 - tb_1,
    bhs_ch = bhs_0 - bhs_1,
    si_ch_r = qnorm((rank(si_ch, na.last = "keep", ties.method = "random") - 0.5) / length(si_ch)), # nolint
  )

lpa_mi_df_mids <-
  imp_t0t1_5 %>%
  complete(., action = "long", include = TRUE) %>%
  group_by(ID) %>%
  mutate(
    si_ch = si_0 - si_1,
    pb_ch = pb_0 - pb_1,
    tb_ch = tb_0 - tb_1,
    bhs_ch = bhs_0 - bhs_1,
    si_ch_r = qnorm((rank(si_ch, na.last = "keep", ties.method = "random") - 0.5) / length(si_ch)), # nolint
  ) %>%
  ungroup() %>%
  as.mids()


lm_mi <- with(lpa_mi_df_mids, lm(si_ch_r ~ pb_ch * tb_ch * bhs_ch))
summary(pool(lm_mi))
tidy(pool(lm_mi), conf.int = T)


lm_ch_class <-
  list(
    full = lm(
      si_ch ~ pb_ch * tb_ch * bhs_ch,
      data = lpa_mi_df,
    ),
    c1 = lm(
      si_ch ~ pb_ch * tb_ch * bhs_ch,
      data = lpa_mi_df[lpa_mi_df$c %in% 1, ],
    ),
    c2 = lm(
      si_ch ~ pb_ch * tb_ch * bhs_ch,
      data = lpa_mi_df[lpa_mi_df$c %in% 2, ],
    ),
    c3 = lm(
      si_ch ~ pb_ch * tb_ch * bhs_ch,
      data = lpa_mi_df[lpa_mi_df$c %in% 3, ],
    )
  )

lm_ch_r_class <-
  list(
    full = lm(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch,
      data = lpa_mi_df,
    ),
    c1 = lm(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch,
      data = lpa_mi_df[lpa_mi_df$c %in% 1, ],
    ),
    c2 = lm(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch,
      data = lpa_mi_df[lpa_mi_df$c %in% 2, ],
    ),
    c3 = lm(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch,
      data = lpa_mi_df[lpa_mi_df$c %in% 3, ],
    )
  )


map(lm_ch_class, ~ summary(.x))
map(lm_ch_class, ~ tidy(.x, conf.int = TRUE, conf.method = "bca"))
map(lm_ch_class, ~ plot(.x))
map(lm_ch_r_class, ~ summary(.x))
map(lm_ch_r_class, ~ tidy(.x, conf.int = TRUE, conf.method = "bca"))
map(lm_ch_r_class, ~ plot(.x))

plot(lm_ch_class$full)
plot(lm_ch_r_class$full) # looks much better
## not convinced much difference here, small sample size
plot(lm_ch_class$c1)
plot(lm_ch_r_class$c1)
## not convinced much difference here
plot(lm_ch_class$c2)
plot(lm_ch_r_class$c2)
## not convinced much difference here
plot(lm_ch_class$c3)
plot(lm_ch_r_class$c3)

save.image("./data/proc/mi_rda.RData")
saveRDS(mget(ls()), "./data/proc/mi_df.rds")
