### ============================================================================
### FILE INFO
### ============================================================================

## func.R: Contains all of the functions needed to perform the actual analysis.
## source()'ing this file should have no side effects other than loading up the
## function definitions. This means that you can modify this file and reload it
## without having to go back an repeat steps 1 & 2 which can take a long time to
## run for large data sets.

### ============================================================================
### LIBRARIES
### ============================================================================

## set seed
addTaskCallback(function(...) {set.seed(46810);TRUE})

## load libraries
library(mice)             # use for multiple imputation
library(miceadds)         # for anova analysis on mids obejcts
## library(sjmisc)        # used to combine mids from mice
## library(tidyLPA)       # for comparison between mclust and mplus
library(MplusAutomation)  # for automating mplus
## NOTE any functions used in analysis should be declarted explicitly
library(aj.HelpRs)        # load my helpRs
## these librares test multivariate normality
## also see http://dwoll.de/rexrepos/posts/normality.html
library(MVN)
library(glue)
library(rhdf5)
library(here)
library(janitor)
library(JTRCI)
library(psych)            # for general analysis
library(lmerTest)
library(lme4)
## library(merTools)
library(car)
library(mlogit)
## library(skimr)
library(lubridate)
library(tidyverse)
## library(viridis)       # viridis colour palette
library(hrbrthemes)       # themes for ggplot2
library(patchwork)
## library(ggtext)
## library(ggpubr)        # pre-formatted stats graphs
library(interactions)
library(patchwork)        # combine multiple plots
## library(formattable)   # formatting of table objects
library(knitr)
library(kableExtra)
library(conflicted)
library(broom.mixed)
library(effectsize)
## generalized linear mixed models
## remotes::install_github("glmmTMB/glmmTMB/glmmTMB")
library(glmmTMB)          # non linear mixed models
options(glmmTMB.cores = 4)
library(bbmle)            # model improvement on glmm
library(DHARMa)           # check overdispersion
library(buildmer)
library(AICcmodavg)
## library(MASS)
library(nnet)
## bootstrapping
library(boot)
library(parallel)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("lmer", "lme4")
conflict_prefer("chisq.test", "stats")

### ============================================================================
### GGPLOT THEMES
### ============================================================================

## ggplot <-
##   function(...)
##     ggplot2::ggplot(...) +
      ## scale_color_brewer(palette = "Set1") +
      ## scale_fill_brewer(palette = "Set1")

### ============================================================================
### KABLE EXTRA OPTION
### ============================================================================
## kableExtra::kbl otptions
kbl <-
  function(...) {
    kableExtra::kbl(
      ...,
      format = getOption("knitr.table.format"),
      booktabs = TRUE,
      longtable = TRUE
      ) %>%
      kableExtra::kable_styling(
        latex_options = "scale_down",
        ## latex_options = c("striped", "scale_down", "hold_position", "repeat_header"),
        bootstrap_options = c(
          "striped", "hover", "condensed", "responsive"
        ),
        full_width = FALSE
      )
  }

### ============================================================================
### LOCAL FUNCTIONS
### ============================================================================

## negation of %in% operator
`%!in%` <- Negate(`%in%`)

## function to quickly create tmp files based on time of save
.tmp_csv <-
  function(data, x) {
    readr::write_csv(data,
      paste0(
        "./data/tmp/",
        ## name writer as string i.e., "name"
        x,
        format(Sys.time(), "_%Y-%m-%d_%H.%M.%S"),
        ".csv"
      ),
      na = ""
    )
  }

## function to quickly final rds for later analysis
.final_rds <-
  function(object, x) {
    saveRDS(object,
      paste0(
        "./data/proc/",
        ## name writer as string i.e., "name"
        x,
        format(Sys.time(), "_%Y_%m_%d"),
        ".rds"
      )
    )
  }

## function to collapse duplicates entries
.coalesce_by_column <-
  function(df) {
    return(dplyr::coalesce(!!!as.list(df)))
  }

## NOTE usage, do not run
## df <-
##   df %>%
##   group_by(ID) %>%
##   summarise_all(coalesce_by_column) %>%
##   ungroup()

## calculated pro rated means and percent missing
.pro_mean_miss <-
  function(data, new_name, ...) {
    data <-
      data %>%
      dplyr::mutate(
        ## mean of rows
        UQ(paste(rlang::syms(c(new_name)), "mean", sep = "_")) :=
          round(
            base::rowMeans(across(...), na.rm = TRUE) *
              ## scaled up to base measures
              base::rowSums(!is.na(across(...))),
            digits = 2
          ),
        ## percent missing of rows
        UQ(paste(rlang::syms(c(new_name)), "miss", sep = "_")) :=
          base::rowMeans(is.na(across(...)))
      )
    return(data)
  }


## calculated pro rated means and percent missing
.pro_mean <-
  function(data, new_name, ...) {
    data <-
      data %>%
      dplyr::mutate(
        ## mean of rows
        UQ(paste(rlang::syms(new_name))) :=
          if_else(
            ## if less than or equal to 20% missing
            base::rowMeans(is.na(across(...))) <= 0.2,
            ## calucate prorated mean
            round(
              base::rowMeans(across(...), na.rm = TRUE) *
                ## scaled up to base measures
                (base::rowSums(is.na(across(...))) +
                  base::rowSums(!is.na(across(...)))),
              digits = 2
            ),
            ## else return NA
            NA_real_
          )
      )
    return(data)
  }
## RELIABLE AND SIGNIFICANT CHANGE
## =============================================================================
## reliable change
.rci_jt <-
  function(time_1, time_2) {
    x <-
      ## (time_1 - time_2)/Sdiff
      (time_2 - time_1) /
      ## Sdiff = sqrt(2 * SE^2)
      sqrt(
        2 * (
          ## SE = sd(time_1) * sqrt(1 - test-retest reliability)
          ## SE = sd(time_1) * sqrt(1 - cor(time_1, time_2))
          exp(
            sd(time_1, na.rm = TRUE) * sqrt(1 - cor(time_1, time_2)))
        )
      )
    return(x)
  }

.sigch_jt <-
  function(time_2, var_rc, var_co) {
    x <-
      case_when(
        var_rc <  -1.96 & time_2 <= var_co ~ "rec",
        var_rc >= -1.96 & time_2 <= var_co ~ "imp",
        ## TRUE ~ "noimp_det"
        var_rc <  -1.96 & time_2 > var_co ~ "norec_sigch",
        var_rc >= -1.96 & time_2 > var_co ~ "noimp_det",
      )
    return(x)
  }

## EFA SCREE PLOT
## -----------------------------------------------------------------------------
## extract relavent data & prepare dataframe for plot
## NOTE requires the output of an mplus object e.g.,
## efa_summary <- readModels(here("efa_dir", "mplus_file.out"))
.efa_scree <-
  function(efa_summary, .title = "EFA Scree Plot") {
    efa_title <- .title
    x <-
      list(
        EFA = efa_summary[["gh5"]][["efa"]][["eigenvalues"]],
        Parallel = efa_summary[["gh5"]][["efa"]][["parallel_average"]]
      )
    ## this
    plot_data <- as_data_frame(x)
    plot_data <- cbind(Factor = paste0(1:nrow(plot_data)), plot_data)
    ## then this
    plot_data <-
      plot_data %>%
      mutate(Factor = fct_inorder(Factor))
    ## pivot the dataframe to “long” format
    plot_data_long <-
      plot_data %>%
      pivot_longer(
        EFA:Parallel,             # The columns I'm gathering together
        names_to = "Analysis",    # new column name for existing names
        values_to = "Eigenvalues" # new column name to store values
      )
    ## plot using ggplot
    efa_scree_plot <-
      plot_data_long %>%
      ggplot(aes(
        y = Eigenvalues,
        x = Factor,
        group = Analysis,
        color = Analysis
      )) +
      geom_point() +
      geom_line() +
      labs(title = efa_title)
      theme_minimal()
    return(efa_scree_plot)
  }

## mice inspect missingness
.miss_df <-
  function(df) {
    .md_df <-
      df %>%
      mice::flux()
    .md_df_names <-
      names(df) %>%
      tibble::as_tibble() %>%
      dplyr::rename(var_name = 1)
    .md_summary <-
      dplyr::bind_cols(
        .md_df_names,
        .md_df
      )
    return(.md_summary)
  }

## mice inspect imputation method per var
.mice_meth_df <-
  function(df) {
    .meth_df <-
      df %>%
      as_tibble()
    .meth_df_names <-
      attr(df, "names") %>%
      tibble::as_tibble() %>%
      dplyr::rename(var_name = 1)
    .meth_summary <-
      dplyr::bind_cols(
        .meth_df_names,
        .meth_df
      ) %>%
      rename(old_value = value) %>%
      mutate(
        new_value = old_value,
        ID = row_number()
      )
    return(.meth_summary)
  }

## mice inspect imputation method per var, the oppositie to the .mice_meth_df
## and expects the object to be in that form
.mice_meth_matrix <-
  function(df) {
    values <-
      df %>%
      mutate(across(everything(), ~ replace_na(.x, ""))) %>%
      select(new_value) %>%
      deframe()
    value_names <-
      df %>%
      mutate(across(everything(), ~ replace_na(.x, ""))) %>%
      select(var_name) %>%
      deframe()
    value_names
    names(values) <- value_names
    return(values)
  }

## mice inspect predictor matrix
.mice_pred_df <-
  function(df) {
    .pred_df <-
      df %>%
      as_tibble()
    .pred_df_names <-
      attr(df, "dimnames")[[1]] %>%
      tibble::as_tibble() %>%
      dplyr::rename(var_name = 1)
    .pred_summary <-
      dplyr::bind_cols(
        .pred_df_names,
        .pred_df
      )
    return(.pred_summary)
  }

## return tibble to pred matrix to be used for mice
.mice_pred_matrix <-
  function(df) {
    .pred_matrix <-
      df %>%
      ## drops the "vars" column
      dplyr::select(-1) %>%
      as.matrix()
    attr(.pred_matrix, "dimnames")[[1]] <-
      attr(.pred_matrix, "dimnames")[[2]]
    return(.pred_matrix)
  }

## symmary by grouping var
## https://www.natedayta.com/2018/03/04/split-a-tidyverse-incarnation-of-split/
.split_ <-
  function(data, ..., .drop = TRUE) {
    stopifnot(inherits(data, "data.frame"))
    vars <-
      ensyms(...)
    vars <-
      purrr::map(
        vars,
        function(x) factor(rlang::eval_tidy(x, data), exclude = NULL)
      )
    base::split(data, vars, drop = .drop)
  }

.summary_by <-
  function(data, ..., .drop = TRUE) {
    stopifnot(inherits(data, "data.frame"))
    vars <-
      ensyms(...)
    vars <-
      purrr::map(
        vars,
        function(x) factor(rlang::eval_tidy(x, data), exclude = NULL)
      )
    new_df <-
      base::split(data, vars, drop = .drop)
    summary_by <-
      new_df %>%
      purrr::map(base::summary)
    return(summary_by)
  }

## stackoverflow.com/questions/46983716/does-a-multi-value-purrrpluck-exist
.pluck_multiple <-
  function(x, ...) {
    `[`(x, ...)
  }

## calculate ICC for lmer model object
icc_lmer <-
  function(m) {
    vc <- as.data.frame((VarCorr(m)))
    l <- vc$vcov
    tibble(grp=vc$grp, icc=sapply(l, function(x){x/sum(l)}))
  }


## by Ben Bolker
## stackoverflow.com/questions/25142901/standardized-coefficients-for-lmer-model
stdCoef.merMod <-
  function(object) {
    sdy <- sd(getME(object,"y"))
    sdx <- apply(getME(object,"X"), 2, sd)
    sc <- fixef(object)*sdx/sdy
    se.fixef <- coef(summary(object))[,"Std. Error"]
    se <- se.fixef*sdx/sdy
    return(data.frame(stdcoef=sc, stdse=se))
  }


## ANOVAS
## -----------------------------------------------------------------------------

## aj_anova <-
##   function(var, group, data) {
##     ## run anovas
##     ## glue("{var}_an") <-
##     ## .args <-
##     ##   rlang::enquos(
##     ##     dv = dv, wid = wid, between = between,
##     ##     within = within, covariate = covariate
##     ##   ) %>%
##     ##   select_quo_variables(data)
##     var_an <-
##       car::Anova(lm(var ~ as.factor(group), data = data), type = 3)
##     ## pairwise t-tests
##     ## var_ttest <-
##     ##   pairwise.t.test(
##     ##     data$var,
##     ##     data$group,
##     ##     p.adjust.method = "bonferroni"
##     ##   )
##     ## ## generate boxplots
##     ## var_box <-
##     ##   data %>%
##     ##   ggplot2::ggplot(aes(as.factor(group), var, fill = factor(group))) +
##     ##   geom_boxplot() +
##     ##   viridis::scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
##     ##   geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
##     ##   hrbrthemes::theme_ipsum() +
##     ##   theme(
##     ##     legend.position = "none",
##     ##     plot.title = element_text(size = 11)
##     ##   ) +
##     ##   ggtitle("A boxplot with jitter") +
##     ##   xlab("")
##     ## ## generate violin plots
##     ## var_violin <-
##     ##   data %>%
##     ##   ggplot(aes(factor(group), var, fill = factor(group))) +
##     ##   geom_violin() +
##     ##   viridis::scale_fill_viridis(discrete = TRUE, alpha = 0.6, option = "A") +
##     ##   hrbrthemes::theme_ipsum() +
##     ##   theme(
##     ##     legend.position = "none",
##     ##     plot.title = element_text(size = 11)
##     ##   ) +
##     ##   ggtitle("PB by group") +
##     ##   xlab("")
##     ## create list of objects
##     var_an_summary <- broom::tidy(var_an)
##     ## broom::glance(var_an)
##     ## summary(var_an)
##     return(var_an_summary)

##   }

## PAIRWISE T TEST
## =============================================================================
## https://stackoverflow.com/questions/27544438/how-to-get-df-and-t-values-from-pairwise-t-test
pairwise.t.test.with.t.and.df <-
  function(x, g, p.adjust.method = p.adjust.methods, pool.sd = !paired,
           paired = FALSE, alternative = c("two.sided", "less", "greater"),
           ...) {
    if (paired & pool.sd) {
      stop("pooling of SD is incompatible with paired tests")
    }
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    g <- factor(g)
    p.adjust.method <- match.arg(p.adjust.method)
    alternative <- match.arg(alternative)
    if (pool.sd) {
      METHOD <- "t tests with pooled SD"
      xbar <- tapply(x, g, mean, na.rm = TRUE)
      s <- tapply(x, g, sd, na.rm = TRUE)
      n <- tapply(!is.na(x), g, sum)
      degf <- n - 1
      total.degf <- sum(degf)
      pooled.sd <- sqrt(sum(s^2 * degf) / total.degf)
      compare.levels <- function(i, j) {
        dif <- xbar[i] - xbar[j]
        se.dif <- pooled.sd * sqrt(1 / n[i] + 1 / n[j])
        t.val <- dif / se.dif
        if (alternative == "two.sided") {
          2 * pt(-abs(t.val), total.degf)
        } else {
          pt(t.val, total.degf, lower.tail = (alternative ==
            "less"))
        }
      }
      compare.levels.t <- function(i, j) {
        dif <- xbar[i] - xbar[j]
        se.dif <- pooled.sd * sqrt(1 / n[i] + 1 / n[j])
        t.val <- dif / se.dif
        t.val
      }
    } else {
      METHOD <- if (paired) {
        "paired t tests"
      } else {
        "t tests with non-pooled SD"
      }
      compare.levels <- function(i, j) {
        xi <- x[as.integer(g) == i]
        xj <- x[as.integer(g) == j]
        t.test(xi, xj,
          paired = paired, alternative = alternative,
          ...
        )$p.value
      }
      compare.levels.t <- function(i, j) {
        xi <- x[as.integer(g) == i]
        xj <- x[as.integer(g) == j]
        t.test(xi, xj,
          paired = paired, alternative = alternative,
          ...
        )$statistic
      }
      compare.levels.df <- function(i, j) {
        xi <- x[as.integer(g) == i]
        xj <- x[as.integer(g) == j]
        t.test(xi, xj,
          paired = paired, alternative = alternative,
          ...
        )$parameter
      }
    }
    PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
    TVAL <- pairwise.table.t(compare.levels.t, levels(g), p.adjust.method)
    if (pool.sd) {
      DF <- total.degf
    } else {
      DF <- pairwise.table.t(compare.levels.df, levels(g), p.adjust.method)
    }
    ans <- list(
      method = METHOD, data.name = DNAME, p.value = PVAL,
      p.adjust.method = p.adjust.method, t.value = TVAL, dfs = DF
    )
    class(ans) <- "pairwise.htest"
    ans
  }
pairwise.table.t <- function(compare.levels.t, level.names, p.adjust.method) {
  ix <- setNames(seq_along(level.names), level.names)
  pp <- outer(ix[-1L], ix[-length(ix)], function(ivec, jvec) {
    sapply(
      seq_along(ivec),
      function(k) {
        i <- ivec[k]
        j <- jvec[k]
        if (i > j) {
          compare.levels.t(i, j)
        } else {
          NA
        }
      }
    )
  })
  pp[lower.tri(pp, TRUE)] <- pp[lower.tri(pp, TRUE)]
  pp
}
### ======================================================================
### END SCRIPT: func.R
### ============================================================================
