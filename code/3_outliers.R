### ============================================================================
### SOURCE SCRIPTS
### ============================================================================
## rm(list = ls(all = TRUE))                                       # nolint
## load("./data/proc/lpa_rda.RData")                               # nolint
## rm(list = setdiff(ls(), c("glm_ch_r_diag", "lpa_t0t1_outc_w"))) # nolint
load("./data/proc/lpa_ol_rda.RData")
source("./code/func.R")
source("./code/load.R")

workspace_size()

### ============================================================================
### SCRATCH OUTLIERS TESTS NOT NEEDED DO NOT KEEP: MOVE TO BAK
### ============================================================================
## OUTLIERS SENSITIVITY ANALYSIS
## =============================================================================
## OUTLIERS TEST 1
## -----------------------------------------------------------------------------
## no significant outliers as tested by DHARMa package
glm_chrol_1_test <-
  map(
    glm_ch_r_diag,
    ~ testOutliers(.x, type = c("bootstrap"), nBoot = 1e3)
  )

map(glm_ch_r_diag, ~ testResiduals(.x))

glm_chrol_1 <-
  map(
    glm_ch_r_diag,
    ## NOTE standard DHARMa outlier test, below more stringent
    ~ outliers(.x,
      ## lowerQuantile = 0.025,
      ## upperQuantile = 0.975
    )
  )

glm_ch_r_ol_1 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_1, ~ summary(.x))
map(glm_ch_r_class, ~ summary(.x))

glm_ch_r_ol_1_dx <-
  map(
    glm_ch_r_ol_1,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_1_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_1,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_1,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_1,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )

## OUTLIER TEST 2
## -----------------------------------------------------------------------------
glm_chrol_2 <-
  map(
    glm_ch_r_ol_1_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_2 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_2, ~ summary(.x))

glm_ch_r_ol_2_dx <-
  map(
    glm_ch_r_ol_2,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_2_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_2,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_2,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_2,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )

## OUTLIERS TEST 3
## -----------------------------------------------------------------------------
glm_chrol_3 <-
  map(
    glm_ch_r_ol_2_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_3 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_3, ~ summary(.x))

glm_ch_r_ol_3_dx <-
  map(
    glm_ch_r_ol_3,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_3_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_3,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_3,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_3,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )

## OUTLIERS TEST 4
## -----------------------------------------------------------------------------
glm_chrol_4 <-
  map(
    glm_ch_r_ol_3_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_4 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_4, ~ summary(.x))

glm_ch_r_ol_4_dx <-
  map(
    glm_ch_r_ol_4,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_4_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_4,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_4,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_4,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )

## OUTLIERS TEST 5
## -----------------------------------------------------------------------------
glm_chrol_5 <-
  map(
    glm_ch_r_ol_4_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_5 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_5, ~ summary(.x))

glm_ch_r_ol_5_dx <-
  map(
    glm_ch_r_ol_5,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_5_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_5,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_5,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_5,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )

## OUTLIERS TEST 6
## -----------------------------------------------------------------------------
glm_chrol_6 <-
  map(
    glm_ch_r_ol_5_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_6 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_6, ~ summary(.x))

glm_ch_r_ol_6_dx <-
  map(
    glm_ch_r_ol_6,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_6_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_6,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_6,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_6,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )

## OUTLIERS TEST 7
## -----------------------------------------------------------------------------
glm_chrol_7 <-
  map(
    glm_ch_r_ol_6_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_7 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_7, ~ summary(.x))

glm_ch_r_ol_7_dx <-
  map(
    glm_ch_r_ol_7,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_7_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_7,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_7,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_7,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )

## OUTLIERS TEST 8
## -----------------------------------------------------------------------------
glm_chrol_8 <-
  map(
    glm_ch_r_ol_7_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_8 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_8, ~ summary(.x))

glm_ch_r_ol_8_dx <-
  map(
    glm_ch_r_ol_8,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_8_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_8,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_8,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_8,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 9
## -----------------------------------------------------------------------------
glm_chrol_9 <-
  map(
    glm_ch_r_ol_8_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_9 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_9, ~ summary(.x))

glm_ch_r_ol_9_dx <-
  map(
    glm_ch_r_ol_9,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_9_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_9,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_9,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_9,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )

## OUTLIERS TEST 10
## -----------------------------------------------------------------------------
glm_chrol_10 <-
  map(
    glm_ch_r_ol_9_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_10 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_10, ~ summary(.x))

glm_ch_r_ol_10_dx <-
  map(
    glm_ch_r_ol_10,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_10_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_10,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_10,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_10,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )

## OUTLIERS TEST 11
## -----------------------------------------------------------------------------
glm_chrol_11 <-
  map(
    glm_ch_r_ol_10_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_11 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_11, ~ summary(.x))

glm_ch_r_ol_11_dx <-
  map(
    glm_ch_r_ol_11,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_11_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_11,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_11,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_11,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )

## OUTLIERS TEST 12
## -----------------------------------------------------------------------------
glm_chrol_12 <-
  map(
    glm_ch_r_ol_11_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_12 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_12, ~ summary(.x))

glm_ch_r_ol_12_dx <-
  map(
    glm_ch_r_ol_12,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_12_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_12,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_12,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_12,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 13
## -----------------------------------------------------------------------------
glm_chrol_13 <-
  map(
    glm_ch_r_ol_12_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_13 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_13, ~ summary(.x))

glm_ch_r_ol_13_dx <-
  map(
    glm_ch_r_ol_13,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_13_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_13,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_13,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_13,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 14
## -----------------------------------------------------------------------------
glm_chrol_14 <-
  map(
    glm_ch_r_ol_13_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_14 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_14$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_14, ~ summary(.x))

glm_ch_r_ol_14_dx <-
  map(
    glm_ch_r_ol_14,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_14_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_14,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_14,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_14,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 15
## -----------------------------------------------------------------------------
glm_chrol_15 <-
  map(
    glm_ch_r_ol_14_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_15 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_15, ~ summary(.x))

glm_ch_r_ol_15_dx <-
  map(
    glm_ch_r_ol_15,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_15_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_15,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_15,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_15,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )

## OUTLIERS TEST 16
## -----------------------------------------------------------------------------
glm_chrol_16 <-
  map(
    glm_ch_r_ol_15_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_16 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_16, ~ summary(.x))

glm_ch_r_ol_16_dx <-
  map(
    glm_ch_r_ol_16,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_16_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_16,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_16,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_16,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 17
## -----------------------------------------------------------------------------
glm_chrol_17 <-
  map(
    glm_ch_r_ol_16_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_17 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_17, ~ summary(.x))

glm_ch_r_ol_17_dx <-
  map(
    glm_ch_r_ol_17,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_17_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_17,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_17,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_17,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 18
## -----------------------------------------------------------------------------
glm_chrol_18 <-
  map(
    glm_ch_r_ol_17_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_18 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_18, ~ summary(.x))

glm_ch_r_ol_18_dx <-
  map(
    glm_ch_r_ol_18,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_18_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_18,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_18,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_18,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 19
## -----------------------------------------------------------------------------
glm_chrol_19 <-
  map(
    glm_ch_r_ol_18_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_19 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_19, ~ summary(.x))

glm_ch_r_ol_19_dx <-
  map(
    glm_ch_r_ol_19,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_19_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_19,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_19,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_19,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 20
## -----------------------------------------------------------------------------
glm_chrol_20 <-
  map(
    glm_ch_r_ol_19_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_20 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_20, ~ summary(.x))

glm_ch_r_ol_20_dx <-
  map(
    glm_ch_r_ol_20,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_20_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_20,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_20,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_20,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 21
## -----------------------------------------------------------------------------
glm_chrol_21 <-
  map(
    glm_ch_r_ol_20_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_21 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_21, ~ summary(.x))

glm_ch_r_ol_21_dx <-
  map(
    glm_ch_r_ol_21,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_21_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_21,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_21,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_21,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 22
## -----------------------------------------------------------------------------
glm_chrol_22 <-
  map(
    glm_ch_r_ol_21_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_22 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_22, ~ summary(.x))

glm_ch_r_ol_22_dx <-
  map(
    glm_ch_r_ol_22,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_22_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_22,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_22,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_22,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 23
## -----------------------------------------------------------------------------
glm_chrol_23 <-
  map(
    glm_ch_r_ol_22_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_23 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_23, ~ summary(.x))

glm_ch_r_ol_23_dx <-
  map(
    glm_ch_r_ol_23,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_23_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_23,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_23,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_23,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 24
## -----------------------------------------------------------------------------
glm_chrol_24 <-
  map(
    glm_ch_r_ol_23_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_24 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_24, ~ summary(.x))

glm_ch_r_ol_24_dx <-
  map(
    glm_ch_r_ol_24,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_24_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_24,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_24,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_24,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 25
## -----------------------------------------------------------------------------
glm_chrol_25 <-
  map(
    glm_ch_r_ol_24_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_25 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_25, ~ summary(.x))

glm_ch_r_ol_25_dx <-
  map(
    glm_ch_r_ol_25,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_25_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_25,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_25,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_25,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 26
## -----------------------------------------------------------------------------
glm_chrol_26 <-
  map(
    glm_ch_r_ol_25_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_26 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_26, ~ summary(.x))

glm_ch_r_ol_26_dx <-
  map(
    glm_ch_r_ol_26,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_26_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_26,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_26,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_26,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 27
## -----------------------------------------------------------------------------
glm_chrol_27 <-
  map(
    glm_ch_r_ol_26_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_27 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_27, ~ summary(.x))

glm_ch_r_ol_27_dx <-
  map(
    glm_ch_r_ol_27,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_27_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_27,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_27,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_27,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 28
## -----------------------------------------------------------------------------
glm_chrol_28 <-
  map(
    glm_ch_r_ol_27_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_28 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_28, ~ summary(.x))

glm_ch_r_ol_28_dx <-
  map(
    glm_ch_r_ol_28,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_28_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_28,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_28,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_28,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 29
## -----------------------------------------------------------------------------
glm_chrol_29 <-
  map(
    glm_ch_r_ol_28_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_29 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_29, ~ summary(.x))

glm_ch_r_ol_29_dx <-
  map(
    glm_ch_r_ol_29,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_29_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_29,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_29,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_29,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 30
## -----------------------------------------------------------------------------
glm_chrol_30 <-
  map(
    glm_ch_r_ol_29_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_30 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_30, ~ summary(.x))

glm_ch_r_ol_30_dx <-
  map(
    glm_ch_r_ol_30,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_30_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_30,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_30,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_30,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 31
## -----------------------------------------------------------------------------
glm_chrol_31 <-
  map(
    glm_ch_r_ol_30_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_31 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_31, ~ summary(.x))

glm_ch_r_ol_31_dx <-
  map(
    glm_ch_r_ol_31,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_31_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_31,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_31,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_31,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 32
## -----------------------------------------------------------------------------
glm_chrol_32 <-
  map(
    glm_ch_r_ol_31_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_32 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_32, ~ summary(.x))

glm_ch_r_ol_32_dx <-
  map(
    glm_ch_r_ol_32,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_32_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_32,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_32,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_32,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 33
## -----------------------------------------------------------------------------
glm_chrol_33 <-
  map(
    glm_ch_r_ol_32_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_33 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_33, ~ summary(.x))

glm_ch_r_ol_33_dx <-
  map(
    glm_ch_r_ol_33,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_33_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_33,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_33,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_33,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 34
## -----------------------------------------------------------------------------
glm_chrol_34 <-
  map(
    glm_ch_r_ol_33_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_34 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_34, ~ summary(.x))

glm_ch_r_ol_34_dx <-
  map(
    glm_ch_r_ol_34,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_34_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_34,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_34,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_34,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 35
## -----------------------------------------------------------------------------
glm_chrol_35 <-
  map(
    glm_ch_r_ol_34_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_35 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_35, ~ summary(.x))

glm_ch_r_ol_35_dx <-
  map(
    glm_ch_r_ol_35,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_35_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_35,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_35,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_35,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 36
## -----------------------------------------------------------------------------
glm_chrol_36 <-
  map(
    glm_ch_r_ol_35_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_36 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ][-glm_chrol_36$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_36, ~ summary(.x))

glm_ch_r_ol_36_dx <-
  map(
    glm_ch_r_ol_36,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_36_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_36,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_36,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_36,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 37
## -----------------------------------------------------------------------------
glm_chrol_37 <-
  map(
    glm_ch_r_ol_36_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_37 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ][-glm_chrol_36$full, ][-glm_chrol_37$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_37, ~ summary(.x))

glm_ch_r_ol_37_dx <-
  map(
    glm_ch_r_ol_37,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_37_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_37,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_37,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_37,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 38
## -----------------------------------------------------------------------------
glm_chrol_38 <-
  map(
    glm_ch_r_ol_37_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_38 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ][-glm_chrol_36$full, ][-glm_chrol_37$full, ][-glm_chrol_38$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_38, ~ summary(.x))

glm_ch_r_ol_38_dx <-
  map(
    glm_ch_r_ol_38,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_38_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_38,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_38,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_38,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 39
## -----------------------------------------------------------------------------
glm_chrol_39 <-
  map(
    glm_ch_r_ol_38_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_39 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ][-glm_chrol_36$full, ][-glm_chrol_37$full, ][-glm_chrol_38$full, ][-glm_chrol_39$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_39, ~ summary(.x))

glm_ch_r_ol_39_dx <-
  map(
    glm_ch_r_ol_39,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_39_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_39,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_39,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_39,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 40
## -----------------------------------------------------------------------------
glm_chrol_40 <-
  map(
    glm_ch_r_ol_39_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_40 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ][-glm_chrol_36$full, ][-glm_chrol_37$full, ][-glm_chrol_38$full, ][-glm_chrol_39$full, ][-glm_chrol_40$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_40, ~ summary(.x))

glm_ch_r_ol_40_dx <-
  map(
    glm_ch_r_ol_40,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_40_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_40,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_40,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_40,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 41
## -----------------------------------------------------------------------------
glm_chrol_41 <-
  map(
    glm_ch_r_ol_40_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_41 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ][-glm_chrol_36$full, ][-glm_chrol_37$full, ][-glm_chrol_38$full, ][-glm_chrol_39$full, ][-glm_chrol_40$full, ][-glm_chrol_41$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_41, ~ summary(.x))

glm_ch_r_ol_41_dx <-
  map(
    glm_ch_r_ol_41,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_41_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_41,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_41,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_41,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 42
## -----------------------------------------------------------------------------
glm_chrol_42 <-
  map(
    glm_ch_r_ol_41_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_42 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ][-glm_chrol_36$full, ][-glm_chrol_37$full, ][-glm_chrol_38$full, ][-glm_chrol_39$full, ][-glm_chrol_40$full, ][-glm_chrol_41$full, ][-glm_chrol_42$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_42, ~ summary(.x))

glm_ch_r_ol_42_dx <-
  map(
    glm_ch_r_ol_42,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_42_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_42,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_42,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_42,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 43
## -----------------------------------------------------------------------------
glm_chrol_43 <-
  map(
    glm_ch_r_ol_42_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_43 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ][-glm_chrol_36$full, ][-glm_chrol_37$full, ][-glm_chrol_38$full, ][-glm_chrol_39$full, ][-glm_chrol_40$full, ][-glm_chrol_41$full, ][-glm_chrol_42$full, ][-glm_chrol_43$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_43, ~ summary(.x))

glm_ch_r_ol_43_dx <-
  map(
    glm_ch_r_ol_43,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_43_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_43,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_43,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_43,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 44
## -----------------------------------------------------------------------------
glm_chrol_44 <-
  map(
    glm_ch_r_ol_43_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_44 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ][-glm_chrol_36$full, ][-glm_chrol_37$full, ][-glm_chrol_38$full, ][-glm_chrol_39$full, ][-glm_chrol_40$full, ][-glm_chrol_41$full, ][-glm_chrol_42$full, ][-glm_chrol_43$full, ][-glm_chrol_44$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_44, ~ summary(.x))

glm_ch_r_ol_44_dx <-
  map(
    glm_ch_r_ol_44,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_44_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_44,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_44,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_44,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 45
## -----------------------------------------------------------------------------
glm_chrol_45 <-
  map(
    glm_ch_r_ol_44_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_45 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ][-glm_chrol_36$full, ][-glm_chrol_37$full, ][-glm_chrol_38$full, ][-glm_chrol_39$full, ][-glm_chrol_40$full, ][-glm_chrol_41$full, ][-glm_chrol_42$full, ][-glm_chrol_43$full, ][-glm_chrol_44$full, ][-glm_chrol_45$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_45, ~ summary(.x))

glm_ch_r_ol_45_dx <-
  map(
    glm_ch_r_ol_45,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_45_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_45,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_45,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_45,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 46
## -----------------------------------------------------------------------------
glm_chrol_46 <-
  map(
    glm_ch_r_ol_45_dx,
    ~ outliers(.x)
  )

glm_ol_46_df <-
  lpa_t0t1_outc_w[-glm_chrol_1$full, ][-glm_chrol_2$full, ][-glm_chrol_3$full, ][-glm_chrol_4$full, ][-glm_chrol_5$full, ][-glm_chrol_6$full, ][-glm_chrol_7$full, ][-glm_chrol_8$full, ][-glm_chrol_9$full, ][-glm_chrol_10$full, ][-glm_chrol_11$full, ][-glm_chrol_12$full, ][-glm_chrol_13$full, ][-glm_chrol_14$full, ][-glm_chrol_15$full, ][-glm_chrol_16$full, ][-glm_chrol_17$full, ][-glm_chrol_18$full, ][-glm_chrol_19$full, ][-glm_chrol_20$full, ][-glm_chrol_21$full, ][-glm_chrol_22$full, ][-glm_chrol_23$full, ][-glm_chrol_24$full, ][-glm_chrol_25$full, ][-glm_chrol_26$full, ][-glm_chrol_27$full, ][-glm_chrol_28$full, ][-glm_chrol_29$full, ][-glm_chrol_30$full, ][-glm_chrol_31$full, ][-glm_chrol_32$full, ][-glm_chrol_33$full, ][-glm_chrol_34$full, ][-glm_chrol_35$full, ][-glm_chrol_36$full, ][-glm_chrol_37$full, ][-glm_chrol_38$full, ][-glm_chrol_39$full, ][-glm_chrol_40$full, ][-glm_chrol_41$full, ][-glm_chrol_42$full, ][-glm_chrol_43$full, ][-glm_chrol_44$full, ][-glm_chrol_45$full, ] # nolint

glm_ch_r_ol_46 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = glm_ol_46_df[-glm_chrol_46$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_46, ~ summary(.x))

glm_ch_r_ol_46_dx <-
  map(
    glm_ch_r_ol_46,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_46_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_46,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_46,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_46,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 47
## -----------------------------------------------------------------------------
glm_chrol_47 <-
  map(
    glm_ch_r_ol_46_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_47 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = glm_ol_46_df[-glm_chrol_46$full, ][-glm_chrol_47$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_47, ~ summary(.x))

glm_ch_r_ol_47_dx <-
  map(
    glm_ch_r_ol_47,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_47_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_47,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_47,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_47,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 48
## -----------------------------------------------------------------------------
glm_chrol_48 <-
  map(
    glm_ch_r_ol_47_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_48 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = glm_ol_46_df[-glm_chrol_46$full, ][-glm_chrol_47$full, ][-glm_chrol_48$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_48, ~ summary(.x))

glm_ch_r_ol_48_dx <-
  map(
    glm_ch_r_ol_48,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_48_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_48,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_48,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_48,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 49
## -----------------------------------------------------------------------------
glm_chrol_49 <-
  map(
    glm_ch_r_ol_48_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_49 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = glm_ol_46_df[-glm_chrol_46$full, ][-glm_chrol_47$full, ][-glm_chrol_48$full, ][-glm_chrol_49$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_49, ~ summary(.x))

glm_ch_r_ol_49_dx <-
  map(
    glm_ch_r_ol_49,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_49_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_49,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_49,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_49,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 50
## -----------------------------------------------------------------------------
glm_chrol_50 <-
  map(
    glm_ch_r_ol_49_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_50 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = glm_ol_46_df[-glm_chrol_46$full, ][-glm_chrol_47$full, ][-glm_chrol_48$full, ][-glm_chrol_49$full, ][-glm_chrol_50$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_50, ~ summary(.x))

glm_ch_r_ol_50_dx <-
  map(
    glm_ch_r_ol_50,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_50_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_50,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_50,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_50,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 51
## -----------------------------------------------------------------------------
glm_chrol_51 <-
  map(
    glm_ch_r_ol_50_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_51 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = glm_ol_46_df[-glm_chrol_46$full, ][-glm_chrol_47$full, ][-glm_chrol_48$full, ][-glm_chrol_49$full, ][-glm_chrol_50$full, ][-glm_chrol_51$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_51, ~ summary(.x))

glm_ch_r_ol_51_dx <-
  map(
    glm_ch_r_ol_51,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_51_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_51,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_51,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_51,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 52
## -----------------------------------------------------------------------------
glm_chrol_52 <-
  map(
    glm_ch_r_ol_51_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_52 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = glm_ol_46_df[-glm_chrol_46$full, ][-glm_chrol_47$full, ][-glm_chrol_48$full, ][-glm_chrol_49$full, ][-glm_chrol_50$full, ][-glm_chrol_51$full, ][-glm_chrol_52$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_52, ~ summary(.x))

glm_ch_r_ol_52_dx <-
  map(
    glm_ch_r_ol_52,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_52_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_52,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_52,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_52,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 53
## -----------------------------------------------------------------------------
glm_chrol_53 <-
  map(
    glm_ch_r_ol_52_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_53 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = glm_ol_46_df[-glm_chrol_46$full, ][-glm_chrol_47$full, ][-glm_chrol_48$full, ][-glm_chrol_49$full, ][-glm_chrol_50$full, ][-glm_chrol_51$full, ][-glm_chrol_52$full, ][-glm_chrol_53$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_53, ~ summary(.x))

glm_ch_r_ol_53_dx <-
  map(
    glm_ch_r_ol_53,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_53_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_53,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_53,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_53,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 54
## -----------------------------------------------------------------------------
glm_chrol_54 <-
  map(
    glm_ch_r_ol_53_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_54 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = glm_ol_46_df[-glm_chrol_46$full, ][-glm_chrol_47$full, ][-glm_chrol_48$full, ][-glm_chrol_49$full, ][-glm_chrol_50$full, ][-glm_chrol_51$full, ][-glm_chrol_52$full, ][-glm_chrol_53$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_54, ~ summary(.x))

glm_ch_r_ol_54_dx <-
  map(
    glm_ch_r_ol_54,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_54_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_54,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_54,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_54,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 55
## -----------------------------------------------------------------------------
glm_chrol_55 <-
  map(
    glm_ch_r_ol_54_dx,
    ~ outliers(.x)
  )

glm_ch_r_ol_55 <-
  list(
    full = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = glm_ol_46_df[-glm_chrol_46$full, ][-glm_chrol_47$full, ][-glm_chrol_48$full, ][-glm_chrol_49$full, ][-glm_chrol_50$full, ][-glm_chrol_51$full, ][-glm_chrol_52$full, ][-glm_chrol_53$full, ][-glm_chrol_54$full, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c1 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 1, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c2 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      ## no outliers
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 2, ][-glm_chrol_1$c2, ],
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    ),
    c3 = glmmTMB(
      si_ch_r ~ pb_ch * tb_ch * bhs_ch + (1 | ID),
      na.action = na.omit,
      data = lpa_t0t1_outc_w[lpa_t0t1_outc_w$c %in% 3, ][-glm_chrol_1$c3, ][-glm_chrol_2$c3, ][-glm_chrol_3$c3, ][-glm_chrol_4$c3, ][-glm_chrol_5$c3, ][-glm_chrol_6$c3, ][-glm_chrol_7$c3, ][-glm_chrol_8$c3, ][-glm_chrol_9$c3, ][-glm_chrol_10$c3, ][-glm_chrol_11$c3, ][-glm_chrol_12$c3, ][-glm_chrol_13$c3, ][-glm_chrol_13$c3, ][-glm_chrol_15$c3, ][-glm_chrol_16$c3, ][-glm_chrol_17$c3, ][-glm_chrol_18$c3, ][-glm_chrol_19$c3, ][-glm_chrol_20$c3, ][-glm_chrol_21$c3, ][-glm_chrol_22$c3, ][-glm_chrol_23$c3, ][-glm_chrol_24$c3, ][-glm_chrol_25$c3, ][-glm_chrol_26$c3, ][-glm_chrol_27$c3, ][-glm_chrol_28$c3, ][-glm_chrol_29$c3, ][-glm_chrol_30$c3, ][-glm_chrol_31$c3, ][-glm_chrol_32$c3, ][-glm_chrol_33$c3, ][-glm_chrol_34$c3, ][-glm_chrol_35$c3, ][-glm_chrol_36$c3, ][-glm_chrol_37$c3, ][-glm_chrol_38$c3, ][-glm_chrol_39$c3, ][-glm_chrol_40$c3, ], # nolint
      family = gaussian,
      ziformula = ~0,
      REML = FALSE,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 30000, eval.max = 40000),
        profile = TRUE,
      )
    )
  )

map(glm_ch_r_ol_55, ~ summary(.x))

glm_ch_r_ol_55_dx <-
  map(
    glm_ch_r_ol_55,
    ~ simulateResiduals(.x)
  )

glm_ch_r_ol_55_dx_plots <-
  list(
    qq = map(
      glm_ch_r_ol_55,
      ~ wrap_elements(plot = ~ plotQQunif(.x), clip = TRUE)
    ),
    resid = map(
      glm_ch_r_ol_55,
      ~ wrap_elements(plot = ~ plotResiduals(.x), clip = TRUE)
    ),
    hist = map(
      glm_ch_r_ol_55,
      ~ wrap_elements(plot = ~ hist(.x), clip = TRUE)
    )
  ) %>%
  map(
    ~ .x %>%
      wrap_plots() +
      plot_annotation(tag_levels = list(c("FULL", "C1", "C2", "C3")))
  )


## OUTLIERS TEST 56
## -----------------------------------------------------------------------------
glm_chrol_56 <-
  map(
    glm_ch_r_ol_55_dx,
    ~ outliers(.x)
  )

## TODO FOR TOMORROW
## -----------------------------------------------------------------------------
## complete the scatterplots as suggested by Brad to look for outliers

lpa_t0t1_outc_w %>%
  select(contains("_ch"), -c("si_ch_r", "dep_ch"), c, ID) %>%
  pivot_longer(-c(c, si_ch, ID)) %>%
  ggplot(aes(x = value, y = si_ch)) +
  ## ggplot(aes(x = value, y = si_ch, fill = c, color = c)) +
  facet_wrap(name ~ ., scales = "free") +
  geom_point() +
  geom_smooth() +
  theme_minimal()

lpa_t0t1_outc_w %>%
  select(contains("_ch"), -c("si_ch_r", "dep_ch"), c, ID) %>%
  pivot_longer(-c(c, ID)) %>%
  ggplot(aes(x = value)) +
  ## ggplot(aes(x = value, fill = c, color = c)) +
  facet_wrap(name ~ ., scales = "free") +
  geom_freqpoly() +
  theme_minimal()

lpa_t0t1_outc_w %>%
  select(contains("_ch"), -c("si_ch_r", "dep_ch"), c, ID) %>%
  pivot_longer(-c(c, ID)) %>%
  ## ggplot(aes(x = name, y = value)) +
  ggplot(aes(x = name, y = value, group = c, fill = c)) +
  facet_wrap(name ~ ., scales = "free") +
  geom_boxplot() +
  theme_minimal()

tt <-
  lpa_t0t1_outc_w %>%
  select(ends_with("_ch"), c, -dep_ch, ID) %>%
  mutate(
   si_ch2 = log(si_ch + 36),
   pb_ch2 = log(si_ch + 20),
   tb_ch2 = log(tb_ch + 23),
   bhs_ch2 = log(bhs_ch + 9),
  )

tt %>%
  select(contains("_ch2"), c, ID) %>%
  pivot_longer(-c(c, ID)) %>%
  ggplot(aes(x = name, y = value, group = c, fill = c)) +
  facet_wrap(name ~ ., scales = "free") +
  geom_boxplot() +
  theme_minimal()

summary(lm(si_ch2 ~ pb_ch2 * tb_ch2 * bhs_ch2, data = tt))
summary(lm(si_ch ~ pb_ch * tb_ch * bhs_ch, data = tt))


### ============================================================================
### WRITE IMAGE AND DATA
### ============================================================================
save.image("./data/proc/lpa_ol_rda.RData")
saveRDS(mget(ls()), "./data/proc/lpa_ol_df.rds")

### ============================================================================
### END SCRIPT
### ============================================================================
