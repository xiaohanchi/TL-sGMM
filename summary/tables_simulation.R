rm(list = ls())
library(dplyr)
library(tidyr)
library(openxlsx)
library(purrr)
library(patchwork)
library(readr)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(emmeans)
library(splines)
library(stringr)

select <- dplyr::select

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

## Main simulations ====================
### iid error
source("simulation/iid_error/settings.R")
result.path <- "simulation/iid_error"
source("summary/sim_summary.R")

output_tbl <- list()
for (ii in seq_along(unique(all.config$rho_source))) {
  output_tbl[[ii]] <- as.data.frame(output) %>%
    tibble() %>%
    mutate(Method = output_methods, .after = 1) %>%
    slice(c(uniq_sc[[ii]])) %>%
    mutate(Scenario = rep(settings[[ii]]$Scenario, each = nrow(uniq_sc[[ii]])), .before = 1) %>%
    select(-c(index)) %>%
    left_join(settings[[ii]], by = c("Scenario")) %>%
    relocate(names(settings[[ii]]), .after = Scenario)

  output_tbl[[ii]] <- output_tbl[[ii]] %>%
    mutate(Scenario = case_when(
      Scenario == 1 ~ 1,
      TRUE ~ floor((Scenario + 1) / 3)
    )) %>%
    filter(re_factor != 1)
}

output_iid <- output_tbl

### ar error
source("simulation/ar_error/settings.R")
result.path <- "simulation/ar_error"
source("summary/sim_summary.R")

output_tbl <- list()
for (ii in seq_along(unique(all.config$rho_source))) {
  output_tbl[[ii]] <- as.data.frame(output) %>%
    tibble() %>%
    mutate(Method = output_methods, .after = 1) %>%
    slice(c(uniq_sc[[ii]])) %>%
    mutate(Scenario = rep(settings[[ii]]$Scenario, each = nrow(uniq_sc[[ii]])), .before = 1) %>%
    select(-c(index)) %>%
    left_join(settings[[ii]], by = c("Scenario")) %>%
    relocate(names(settings[[ii]]), .after = Scenario)

  output_tbl[[ii]] <- output_tbl[[ii]] %>%
    mutate(Scenario = (floor((Scenario - 1) / 3) + 1)) %>%
    filter(re_factor != 1)
}

output_ar <- output_tbl

## Tables: Paper ======================================
### Simulation ==========================
sample_size <- 80
#### Table 1 ===========================================
display_col <- c("Scenario", "Method", "fpr", "fnr", "cfpr_3sigma", "cfpr_2sigma", "cfnr_3sigma", "cfnr_2sigma")
table1 <- list()
for (ii in 1:2) {
  beta_shift_type <- c(1, 3)[ii]
  table1_tmp <- list()
  ### iid error
  rho_source <- 0.5
  list_id <- which(sapply(output_iid, function(vv) unique(vv$rho_source)) == rho_source)
  table1_tmp[[1]] <- output_iid[[list_id]] %>%
    filter(beta_delta_idx == beta_shift_type) %>%
    select(all_of(display_col)) %>%
    mutate(
      Method = factor(Method,
        levels = c("onsite_smm", "oracle_smm", "pool_smm", "2stage_smm", "rsmm", "gmm3", "gmm4", "gmm5")
      )
    ) %>%
    mutate(Method = recode(Method,
      "onsite_smm" = "Onsite SMM",
      "oracle_smm" = "Oracle SMM",
      "pool_smm"   = "Pooled SMM",
      "2stage_smm" = "Two-stage SMM",
      "rsmm"       = "Robust SMM",
      "gmm3"       = "Pooled sGMM v0",
      "gmm4"       = "Pooled sGMM",
      "gmm5"       = "TL-sGMM"
    )) %>%
    filter(!(Method %in% c("Onsite SMM", "Oracle SMM", "Pooled sGMM v0")))


  scenario_levels <- sort(unique(table1_tmp[[1]]$Scenario))
  table1_tmp[[1]] <- table1_tmp[[1]] %>%
    mutate(Scenario = factor(Scenario,
      levels = scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
    )) %>%
    mutate(Scenario = recode(
      Scenario,
      !!!setNames(
        c("1", "2a", "2b", "2c", "2d", "3a", "3b"),
        scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
      )
    )) %>%
    arrange(Scenario, Method) %>%
    replace_na(list(fnr = 0, cfnr = 0)) %>%
    mutate(blank1 = NA, .after = 6)

  empty_row <- table1_tmp[[1]][1, ] %>% mutate(across(everything(), ~NA))
  table1_tmp[[1]] <- table1_tmp[[1]] %>%
    group_by(Scenario) %>%
    group_split() %>%
    map_df(~ bind_rows(.x, empty_row)) %>%
    group_by(Scenario) %>%
    mutate(
      Scenario = if_else(row_number() == 1, Scenario, NA_character_)
    ) %>%
    ungroup()

  ### AR error
  rho_source <- 0.5
  list_id <- which(sapply(output_ar, function(vv) unique(vv$rho_source)) == rho_source)
  table1_tmp[[2]] <- output_ar[[list_id]] %>%
    filter(beta_delta_idx == beta_shift_type) %>%
    select(all_of(display_col)) %>%
    mutate(
      Method = factor(Method,
        levels = c("onsite_smm", "oracle_smm", "pool_smm", "2stage_smm", "rsmm", "gmm3", "gmm4", "gmm5")
      )
    ) %>%
    mutate(Method = recode(Method,
      "onsite_smm" = "Onsite SMM",
      "oracle_smm" = "Oracle SMM",
      "pool_smm"   = "Pooled SMM",
      "2stage_smm" = "Two-stage SMM",
      "rsmm"       = "Robust SMM",
      "gmm3"       = "Pooled sGMM v0",
      "gmm4"       = "Pooled sGMM",
      "gmm5"       = "TL-sGMM"
    )) %>%
    filter(!(Method %in% c("Onsite SMM", "Oracle SMM", "Pooled sGMM v0")))


  scenario_levels <- sort(unique(table1_tmp[[2]]$Scenario))
  table1_tmp[[2]] <- table1_tmp[[2]] %>%
    mutate(Scenario = factor(Scenario,
      levels = scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
    )) %>%
    mutate(Scenario = recode(
      Scenario,
      !!!setNames(
        c("1", "2a", "2b", "2c", "2d", "3a", "3b"),
        scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
      )
    )) %>%
    arrange(Scenario, Method) %>%
    replace_na(list(fnr = 0, cfnr = 0)) %>%
    mutate(blank1 = NA, .after = 6)

  empty_row <- table1_tmp[[2]][1, ] %>% mutate(across(everything(), ~NA))
  table1_tmp[[2]] <- table1_tmp[[2]] %>%
    group_by(Scenario) %>%
    group_split() %>%
    map_df(~ bind_rows(.x, empty_row)) %>%
    group_by(Scenario) %>%
    mutate(
      Scenario = if_else(row_number() == 1, Scenario, NA_character_)
    ) %>%
    ungroup()

  table1[[ii]] <- bind_cols(
    table1_tmp[[1]],
    NA,
    (table1_tmp[[2]] %>% select(-c(Scenario, Method)))
  )
}

table1 <- list(
  "consistent" = table1[[1]],
  "shift" = table1[[2]]
)


#### Table 2 ===========================================
display_col <- c("Scenario", "Method", "rbias_point", "rmse_point", "CI_length", "CI_coverage", "rbias_fwc", "rmse_fwc", "CI_length_fwc", "CI_coverage_fwc")

table2 <- list()
for (ii in 1:2) {
  beta_shift_type <- c(1, 3)[ii]
  table2_tmp <- list()
  ### iid error
  list_id <- which(sapply(output_iid, function(vv) unique(vv$rho_source)) == 0.5)
  table2_tmp[[1]] <- output_iid[[list_id]] %>%
    filter(beta_delta_idx == beta_shift_type) %>%
    select(all_of(display_col)) %>%
    mutate(
      Method = factor(Method,
        levels = c("onsite_smm", "oracle_smm", "pool_smm", "2stage_smm", "rsmm", "gmm3", "gmm4", "gmm5")
      )
    ) %>%
    mutate(Method = recode(Method,
      "onsite_smm" = "Onsite SMM",
      "oracle_smm" = "Oracle SMM",
      "pool_smm"   = "Pooled SMM",
      "2stage_smm" = "Two-stage SMM",
      "rsmm"       = "Robust SMM",
      "gmm3"       = "Pooled sGMM v0",
      "gmm4"       = "Pooled sGMM",
      "gmm5"       = "TL-sGMM"
    )) %>%
    filter(!(Method %in% c("Pooled sGMM v0")))

  scenario_levels <- sort(unique(table2_tmp[[1]]$Scenario))
  table2_tmp[[1]] <- table2_tmp[[1]] %>%
    mutate(Scenario = factor(Scenario,
      levels = scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
    )) %>%
    mutate(Scenario = recode(
      Scenario,
      !!!setNames(
        c("1", "2a", "2b", "2c", "2d", "3a", "3b"),
        scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
      )
    )) %>%
    arrange(Scenario, Method) %>%
    mutate(blank2 = NA, .after = 6)

  empty_row <- table2_tmp[[1]][1, ] %>% mutate(across(everything(), ~NA))
  table2_tmp[[1]] <- table2_tmp[[1]] %>%
    group_by(Scenario) %>%
    group_split() %>%
    map_df(~ bind_rows(.x, empty_row)) %>%
    group_by(Scenario) %>%
    mutate(
      Scenario = if_else(row_number() == 1, Scenario, NA_character_)
    ) %>%
    ungroup()

  ### AR error
  list_id <- which(sapply(output_ar, function(vv) unique(vv$rho_source)) == 0.5)
  table2_tmp[[2]] <- output_ar[[list_id]] %>%
    filter(beta_delta_idx == beta_shift_type) %>%
    select(all_of(display_col)) %>%
    mutate(
      Method = factor(Method,
        levels = c("onsite_smm", "oracle_smm", "pool_smm", "2stage_smm", "rsmm", "gmm3", "gmm4", "gmm5")
      )
    ) %>%
    mutate(Method = recode(Method,
      "onsite_smm" = "Onsite SMM",
      "oracle_smm" = "Oracle SMM",
      "pool_smm"   = "Pooled SMM",
      "2stage_smm" = "Two-stage SMM",
      "rsmm"       = "Robust SMM",
      "gmm3"       = "Pooled sGMM v0",
      "gmm4"       = "Pooled sGMM",
      "gmm5"       = "TL-sGMM"
    )) %>%
    filter(!(Method %in% c("Pooled sGMM v0")))

  scenario_levels <- sort(unique(table2_tmp[[2]]$Scenario))
  table2_tmp[[2]] <- table2_tmp[[2]] %>%
    mutate(Scenario = factor(Scenario,
      levels = scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
    )) %>%
    mutate(Scenario = recode(
      Scenario,
      !!!setNames(
        c("1", "2a", "2b", "2c", "2d", "3a", "3b"),
        scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
      )
    )) %>%
    arrange(Scenario, Method) %>%
    mutate(blank2 = NA, .after = 6)

  empty_row <- table2_tmp[[2]][1, ] %>% mutate(across(everything(), ~NA))
  table2_tmp[[2]] <- table2_tmp[[2]] %>%
    group_by(Scenario) %>%
    group_split() %>%
    map_df(~ bind_rows(.x, empty_row)) %>%
    group_by(Scenario) %>%
    mutate(
      Scenario = if_else(row_number() == 1, Scenario, NA_character_)
    ) %>%
    ungroup()


  if (ii == 1) {
    table2 <- list(
      "iid error, consistent" = table2_tmp[[1]],
      "AR error, consistent"  = table2_tmp[[2]]
    )
  } else if (ii == 2) {
    table2 <- c(
      table2,
      list(
        "iid error, shift" = table2_tmp[[1]],
        "AR error, shift"  = table2_tmp[[2]]
      )
    )
  }
}


### Sensitivity: Number of clusters ====================
sensitivity_env <- new.env()
source("simulation/sensitivity/num_clusters/settings.R", local = sensitivity_env)
with(sensitivity_env, {
  result.path <- "simulation/sensitivity/num_clusters"
  all.config <- all.config %>%
    group_by(gmm_ncls) %>%
    mutate(scenario = row_number()) %>%
    ungroup()
  uniq_sc <- lapply(unique(all.config$rho_source), function(rr) {
    sapply(
      unique(all.config$scenario[all.config$rho_source == rr]),
      function(ss) which(all.config$scenario == ss)
    )
  })
  output <- matrix(NA, nrow = nrow(all.config), ncol = 30)
  output_ncls <- numeric(nrow(all.config))
  sc_vars <- c(
    "n_patients", "beta_delta_idx", "rho_source", "re_factor",
    "remote_inflation", "outlier_rate1", "outlier_rate2",
    "outlier_type", "error_type", "uni_modal", "high_freq", "DGM_type"
  )
  all.config <- all.config %>%
    group_by(across(all_of(sc_vars))) %>%
    mutate(sc_id = cur_group_id()) %>%
    ungroup()
  ok_reps <- vector("list", nrow(all.config))
  for (jj in 1:nrow(all.config)) {
    fjj <- paste0(result.path, "/output_", jj, ".csv")
    if (!file.exists(fjj)) next
    djj <- read_csv(fjj, show_col_types = FALSE)
    if (all(grepl("gmm", djj$method))) {
      ok_reps[[jj]] <- djj$replicate[djj$converge_flag == 1]
    } else {
      ok_reps[[jj]] <- djj$replicate
    }
  }
  common_reps <- lapply(split(seq_len(nrow(all.config)), all.config$sc_id), function(idxs) {
    sets <- ok_reps[idxs]
    sets <- sets[!vapply(sets, is.null, logical(1))]
    if (length(sets) == 0) {
      return(NULL)
    }
    Reduce(intersect, sets)
  })

  for (ii in 1:nrow(all.config)) {
    filepath <- paste0(result.path, "/output_", ii, ".csv")
    filepath_muhat <- paste0(result.path, "/muhat_fixed_", ii, ".csv")
    filepath_mutrue <- paste0(result.path, "/mutrue_fixed_", ii, ".csv")
    if (!file.exists(filepath)) next

    tmp_res <- read_csv(filepath, show_col_types = FALSE)
    n_all <- nrow(tmp_res)
    if (all(grepl("gmm", tmp_res$method))) {
      tmp_res <- tmp_res[tmp_res$converge_flag == 1, ]
    }
    reps_keep <- common_reps[[as.character(all.config$sc_id[ii])]]
    if (!is.null(reps_keep)) {
      tmp_res <- tmp_res[tmp_res$replicate %in% reps_keep, ]
    }
    n_conv <- nrow(tmp_res)
    if (n_all < 500) {
      warning(paste0("File output_", ii, " has only ", n_all, " rows (< 500)"))
    }
    fpr <- round(mean(tmp_res$fpr) * 100, 2)
    fnr <- round(mean(tmp_res$fnr) * 100, 2)
    cfpr_3sigma <- round(mean(tmp_res$cfpr_3sigma) * 100, 2)
    cfnr_3sigma <- round(mean(tmp_res$cfnr_3sigma) * 100, 2)
    cfpr_2sigma <- round(mean(tmp_res$cfpr_2sigma) * 100, 2)
    cfnr_2sigma <- round(mean(tmp_res$cfnr_2sigma) * 100, 2)

    error_pts <- round(mean(tmp_res$error_p_num / tmp_res$error_p_denom) * 100, 2)
    cerror_3sigma_pts <- round(mean(tmp_res$cerror_3sigma_p_num / tmp_res$error_p_denom) * 100, 2)
    cerror_2sigma_pts <- round(mean(tmp_res$cerror_2sigma_p_num / tmp_res$error_p_denom) * 100, 2)

    # estimation
    tmp_mu <- list()
    tmp_mu[[1]] <- read_csv(filepath_muhat, show_col_types = FALSE) %>%
      filter(replicate %in% tmp_res$replicate) %>%
      select(-c(replicate, method))
    tmp_mu[[2]] <- read_csv(filepath_mutrue, show_col_types = FALSE) %>%
      filter(replicate %in% tmp_res$replicate) %>%
      select(-c(replicate, method))

    bias_point <- colMeans(tmp_mu[[1]] - tmp_mu[[2]]) %>%
      mean() %>%
      round(2)
    rbias_point <- colMeans((tmp_mu[[1]] - tmp_mu[[2]]) / mean(unlist(distinct(tmp_mu[[2]])))) %>%
      mean() %>%
      {
        round(. * 100, 2)
      }
    rmse_point <- colMeans((tmp_mu[[1]] - tmp_mu[[2]])^2) %>%
      sqrt() %>%
      mean()

    bias_fwc <- mean(tmp_res$fwc_fixed_est - tmp_res$fwc_fixed_true) %>% round(2)
    rbias_fwc <- mean((tmp_res$fwc_fixed_est - tmp_res$fwc_fixed_true) / tmp_res$fwc_fixed_true) %>%
      {
        round(. * 100, 2)
      }
    rmse_fwc <- mean((tmp_res$fwc_fixed_est - tmp_res$fwc_fixed_true)^2) %>% sqrt()

    CI_length <- mean(tmp_res$CI_length) %>% round(3)
    CI_cov <- mean(tmp_res$CI_cov) %>%
      {
        round(. * 100, 2)
      }
    CI_length_fwc <- mean(tmp_res$CI_length_fwc) %>% round(3)
    CI_cov_fwc <- mean(tmp_res$CI_cov_fwc) %>%
      {
        round(. * 100, 2)
      }

    ## in-sample
    rbias_point_insp <- mean(tmp_res$insp_rbias) %>%
      {
        round(. * 100, 2)
      }
    rmse_point_insp <- mean(tmp_res$insp_rmse)
    rbias_fwc_insp <- mean(tmp_res$insp_rbias_fwc) %>%
      {
        round(. * 100, 2)
      }
    rmse_fwc_insp <- mean(tmp_res$insp_rmse_fwc)

    ## out-sample
    rbias_point_outsp <- mean(tmp_res$outsp_rbias) %>%
      {
        round(. * 100, 2)
      }
    rmse_point_outsp <- mean(tmp_res$outsp_rmse)
    rbias_fwc_outsp <- mean(tmp_res$outsp_rbias_fwc) %>%
      {
        round(. * 100, 2)
      }
    rmse_fwc_outsp <- mean(tmp_res$outsp_rmse_fwc)

    output_ncls[ii] <- all.config$gmm_ncls[ii]
    output[ii, ] <- c(
      ii,
      c(
        fpr, fnr, cfpr_3sigma, cfnr_3sigma, cfpr_2sigma, cfnr_2sigma,
        error_pts, cerror_3sigma_pts, cerror_2sigma_pts,
        bias_point, rbias_point, rmse_point,
        bias_fwc, rbias_fwc, rmse_fwc,
        CI_length, CI_cov, CI_length_fwc, CI_cov_fwc,
        rbias_point_insp, rmse_point_insp, rbias_fwc_insp, rmse_fwc_insp,
        rbias_point_outsp, rmse_point_outsp, rbias_fwc_outsp, rmse_fwc_outsp,
        n_conv, n_all
      )
    )
  }
  colnames(output) <- c(
    "index", "fpr", "fnr", "cfpr_3sigma", "cfnr_3sigma", "cfpr_2sigma", "cfnr_2sigma",
    "error_patient", "cerror_3sigma_patient", "cerror_2sigma_patient",
    "bias_point", "rbias_point", "rmse_point", "bias_fwc", "rbias_fwc", "rmse_fwc",
    "CI_length", "CI_coverage", "CI_length_fwc", "CI_coverage_fwc",
    "rbias_point_insp", "rmse_point_insp", "rbias_fwc_insp", "rmse_fwc_insp",
    "rbias_point_outsp", "rmse_point_outsp", "rbias_fwc_outsp", "rmse_fwc_outsp",
    "num_converge", "num_allsimu"
  )
  col.idx <- c("n_patients", "beta_delta_idx", "rho_source", "re_factor", "outlier_rate1", "outlier_rate2", "outlier_type", "error_type", "uni_modal", "high_freq")

  settings <- lapply(seq_along(unique(all.config$rho_source)), function(rr) {
    all.config[uniq_sc[[rr]][1, ], col.idx] %>% mutate(Scenario = 1:n(), .before = 1)
  })

  output_tbl <- list()
  for (ii in seq_along(unique(all.config$rho_source))) {
    output_tbl[[ii]] <- as.data.frame(output) %>%
      tibble() %>%
      mutate(L = output_ncls, .after = 1) %>%
      slice(c(uniq_sc[[ii]])) %>%
      mutate(Scenario = rep(1:nrow(settings[[ii]]), each = nrow(uniq_sc[[ii]])), .before = 1) %>%
      select(-c(index)) %>%
      left_join(settings[[ii]], by = "Scenario") %>%
      relocate(names(settings[[ii]]), .after = Scenario)

    output_tbl[[ii]] <- output_tbl[[ii]] %>%
      mutate(Scenario = (floor((Scenario - 1) / 2) + 1))
  }
})


#### Table S1 ===========================================

display_col <- c("Scenario", "Method", "fpr", "fnr", "cfpr_3sigma", "cfpr_2sigma", "cfnr_3sigma", "cfnr_2sigma")
tableS1 <- list()
tableS1_tmp <- list()
for (ll in 1:2) {
  beta_shift_type <- c(1, 3)[ll]
  tableS1_tmp[[ll]] <- sensitivity_env$output_tbl[[1]] %>%
    filter(beta_delta_idx == beta_shift_type) %>%
    select(all_of(c(setdiff(display_col, "Method"), "L"))) %>%
    mutate(Scenario = factor(Scenario,
      levels = c(1, 4, 5, 2, 3, 7, 6)
    )) %>%
    mutate(Scenario = recode(
      Scenario,
      !!!setNames(
        c("1", "2a", "2b", "2c", "2d", "3a", "3b"),
        c(1, 4, 5, 2, 3, 7, 6)
      )
    )) %>%
    arrange(Scenario, L)

  empty_row <- tableS1_tmp[[ll]][1, ] %>% mutate(across(everything(), ~NA))
  tableS1_tmp[[ll]] <- tableS1_tmp[[ll]] %>%
    group_by(Scenario) %>%
    group_split() %>%
    map_df(~ bind_rows(.x, empty_row)) %>%
    group_by(Scenario) %>%
    mutate(
      Scenario = if_else(row_number() == 1, Scenario, NA_character_)
    ) %>%
    ungroup()

  tableS1_tmp[[ll]] <- tableS1_tmp[[ll]] %>%
    relocate(L, .after = 1) %>%
    mutate(blank1 = NA, .after = 6)
}

tableS1 <- bind_cols(
  tableS1_tmp[[1]],
  NA,
  (tableS1_tmp[[2]] %>% select(-c(Scenario, L)))
)

#### Table S2 ===========================================

display_col <- c("Scenario", "Method", "rbias_point", "rmse_point", "CI_length", "CI_coverage", "rbias_fwc", "rmse_fwc", "CI_length_fwc", "CI_coverage_fwc")

tableS2 <- tableS2_tmp <- list()
for (ll in 1:2) {
  beta_shift_type <- c(1, 3)[ll]
  tableS2_tmp[[ll]] <- sensitivity_env$output_tbl[[1]] %>%
    filter(beta_delta_idx == beta_shift_type) %>%
    select(all_of(c(setdiff(display_col, "Method"), "L"))) %>%
    mutate(Scenario = factor(Scenario,
      levels = c(1, 4, 5, 2, 3, 7, 6)
    )) %>%
    mutate(Scenario = recode(
      Scenario,
      !!!setNames(
        c("1", "2a", "2b", "2c", "2d", "3a", "3b"),
        c(1, 4, 5, 2, 3, 7, 6)
      )
    )) %>%
    arrange(Scenario, L)

  empty_row <- tableS2_tmp[[ll]][1, ] %>% mutate(across(everything(), ~NA))
  tableS2_tmp[[ll]] <- tableS2_tmp[[ll]] %>%
    group_by(Scenario) %>%
    group_split() %>%
    map_df(~ bind_rows(.x, empty_row)) %>%
    group_by(Scenario) %>%
    mutate(
      Scenario = if_else(row_number() == 1, Scenario, NA_character_)
    ) %>%
    ungroup()
  tableS2_tmp[[ll]] <- tableS2_tmp[[ll]] %>%
    relocate(L, .after = 1) %>%
    mutate(blank1 = NA, .after = 6)
}
tableS2 <- bind_rows(tableS2_tmp[[1]], tableS2_tmp[[2]])


### Sensitivity: Number of measures ====================
sensitivity_env <- new.env()
source("simulation/sensitivity/num_measures/settings.R", local = sensitivity_env)
with(sensitivity_env, {
  result.path <- "simulation/sensitivity/num_measures"
  all.config <- all.config %>%
    group_by(model_type, gmm_type) %>%
    mutate(scenario = row_number()) %>%
    ungroup()
  uniq_sc <- lapply(unique(all.config$rho_source), function(rr) {
    sapply(
      unique(all.config$scenario[all.config$rho_source == rr]),
      function(ss) which(all.config$scenario == ss)
    )
  })
  output <- matrix(NA, nrow = nrow(all.config), ncol = 30)
  output_methods <- character(nrow(all.config))
  sc_vars <- c(
    "n_patients", "beta_delta_idx", "rho_source", "re_factor",
    "remote_inflation", "outlier_rate1", "outlier_rate2",
    "outlier_type", "error_type", "uni_modal", "high_freq", "DGM_type"
  )
  all.config <- all.config %>%
    group_by(across(all_of(sc_vars))) %>%
    mutate(sc_id = cur_group_id()) %>%
    ungroup()
  ok_reps <- vector("list", nrow(all.config))
  for (jj in 1:nrow(all.config)) {
    fjj <- paste0(result.path, "/output_", jj, ".csv")
    if (!file.exists(fjj)) next
    djj <- read_csv(fjj, show_col_types = FALSE)
    if (all(grepl("gmm", djj$method))) {
      ok_reps[[jj]] <- djj$replicate[djj$converge_flag == 1]
    } else {
      ok_reps[[jj]] <- djj$replicate
    }
  }
  paper_idx <- which(all.config$model_type %in% 1:5 |
    (all.config$model_type == 6 & all.config$gmm_type %in% c(4, 5)))
  common_reps <- lapply(split(paper_idx, all.config$sc_id[paper_idx]), function(idxs) {
    sets <- ok_reps[idxs]
    sets <- sets[!vapply(sets, is.null, logical(1))]
    if (length(sets) == 0) {
      return(NULL)
    }
    Reduce(intersect, sets)
  })

  for (ii in 1:nrow(all.config)) {
    filepath <- paste0(result.path, "/output_", ii, ".csv")
    filepath_muhat <- paste0(result.path, "/muhat_fixed_", ii, ".csv")
    filepath_mutrue <- paste0(result.path, "/mutrue_fixed_", ii, ".csv")
    if (!file.exists(filepath)) next

    tmp_res <- read_csv(filepath, show_col_types = FALSE)
    n_all <- nrow(tmp_res)
    if (all(grepl("gmm", tmp_res$method))) {
      tmp_res <- tmp_res[tmp_res$converge_flag == 1, ]
    }
    reps_keep <- common_reps[[as.character(all.config$sc_id[ii])]]
    if (!is.null(reps_keep)) {
      tmp_res <- tmp_res[tmp_res$replicate %in% reps_keep, ]
    }
    n_conv <- nrow(tmp_res)
    if (n_all < 500) {
      warning(paste0("File output_", ii, " has only ", n_all, " rows (< 500)"))
    }
    fpr <- round(mean(tmp_res$fpr) * 100, 2)
    fnr <- round(mean(tmp_res$fnr) * 100, 2)
    cfpr_3sigma <- round(mean(tmp_res$cfpr_3sigma) * 100, 2)
    cfnr_3sigma <- round(mean(tmp_res$cfnr_3sigma) * 100, 2)
    cfpr_2sigma <- round(mean(tmp_res$cfpr_2sigma) * 100, 2)
    cfnr_2sigma <- round(mean(tmp_res$cfnr_2sigma) * 100, 2)

    error_pts <- round(mean(tmp_res$error_p_num / tmp_res$error_p_denom) * 100, 2)
    cerror_3sigma_pts <- round(mean(tmp_res$cerror_3sigma_p_num / tmp_res$error_p_denom) * 100, 2)
    cerror_2sigma_pts <- round(mean(tmp_res$cerror_2sigma_p_num / tmp_res$error_p_denom) * 100, 2)

    # estimation
    tmp_mu <- list()
    tmp_mu[[1]] <- read_csv(filepath_muhat, show_col_types = FALSE) %>%
      filter(replicate %in% tmp_res$replicate) %>%
      select(-c(replicate, method))
    tmp_mu[[2]] <- read_csv(filepath_mutrue, show_col_types = FALSE) %>%
      filter(replicate %in% tmp_res$replicate) %>%
      select(-c(replicate, method))

    bias_point <- colMeans(tmp_mu[[1]] - tmp_mu[[2]]) %>%
      mean() %>%
      round(2)
    rbias_point <- colMeans((tmp_mu[[1]] - tmp_mu[[2]]) / mean(unlist(distinct(tmp_mu[[2]])))) %>%
      mean() %>%
      {
        round(. * 100, 2)
      }
    rmse_point <- colMeans((tmp_mu[[1]] - tmp_mu[[2]])^2) %>%
      sqrt() %>%
      mean()

    bias_fwc <- mean(tmp_res$fwc_fixed_est - tmp_res$fwc_fixed_true) %>% round(2)
    rbias_fwc <- mean((tmp_res$fwc_fixed_est - tmp_res$fwc_fixed_true) / tmp_res$fwc_fixed_true) %>%
      {
        round(. * 100, 2)
      }
    rmse_fwc <- mean((tmp_res$fwc_fixed_est - tmp_res$fwc_fixed_true)^2) %>% sqrt()

    CI_length <- mean(tmp_res$CI_length) %>% round(3)
    CI_cov <- mean(tmp_res$CI_cov) %>%
      {
        round(. * 100, 2)
      }
    CI_length_fwc <- mean(tmp_res$CI_length_fwc) %>% round(3)
    CI_cov_fwc <- mean(tmp_res$CI_cov_fwc) %>%
      {
        round(. * 100, 2)
      }

    ## in-sample
    rbias_point_insp <- mean(tmp_res$insp_rbias) %>%
      {
        round(. * 100, 2)
      }
    rmse_point_insp <- mean(tmp_res$insp_rmse)
    rbias_fwc_insp <- mean(tmp_res$insp_rbias_fwc) %>%
      {
        round(. * 100, 2)
      }
    rmse_fwc_insp <- mean(tmp_res$insp_rmse_fwc)

    ## out-sample
    rbias_point_outsp <- mean(tmp_res$outsp_rbias) %>%
      {
        round(. * 100, 2)
      }
    rmse_point_outsp <- mean(tmp_res$outsp_rmse)
    rbias_fwc_outsp <- mean(tmp_res$outsp_rbias_fwc) %>%
      {
        round(. * 100, 2)
      }
    rmse_fwc_outsp <- mean(tmp_res$outsp_rmse_fwc)

    output_methods[ii] <- tmp_res$method[1]
    output[ii, ] <- c(
      ii,
      c(
        fpr, fnr, cfpr_3sigma, cfnr_3sigma, cfpr_2sigma, cfnr_2sigma,
        error_pts, cerror_3sigma_pts, cerror_2sigma_pts,
        bias_point, rbias_point, rmse_point,
        bias_fwc, rbias_fwc, rmse_fwc,
        CI_length, CI_cov, CI_length_fwc, CI_cov_fwc,
        rbias_point_insp, rmse_point_insp, rbias_fwc_insp, rmse_fwc_insp,
        rbias_point_outsp, rmse_point_outsp, rbias_fwc_outsp, rmse_fwc_outsp,
        n_conv, n_all
      )
    )
  }
  colnames(output) <- c(
    "index", "fpr", "fnr", "cfpr_3sigma", "cfnr_3sigma", "cfpr_2sigma", "cfnr_2sigma",
    "error_patient", "cerror_3sigma_patient", "cerror_2sigma_patient",
    "bias_point", "rbias_point", "rmse_point", "bias_fwc", "rbias_fwc", "rmse_fwc",
    "CI_length", "CI_coverage", "CI_length_fwc", "CI_coverage_fwc",
    "rbias_point_insp", "rmse_point_insp", "rbias_fwc_insp", "rmse_fwc_insp",
    "rbias_point_outsp", "rmse_point_outsp", "rbias_fwc_outsp", "rmse_fwc_outsp",
    "num_converge", "num_allsimu"
  )
  col.idx <- c("n_patients", "beta_delta_idx", "rho_source", "re_factor", "outlier_rate1", "outlier_rate2", "outlier_type", "error_type", "uni_modal", "high_freq")

  settings <- lapply(seq_along(unique(all.config$rho_source)), function(rr) {
    all.config[uniq_sc[[rr]][1, ], col.idx] %>% mutate(Scenario = 1:n(), .before = 1)
  })

  output_tbl <- list()
  for (ii in seq_along(unique(all.config$rho_source))) {
    output_tbl[[ii]] <- as.data.frame(output) %>%
      tibble() %>%
      mutate(Method = output_methods, .after = 1) %>%
      slice(c(uniq_sc[[ii]])) %>%
      mutate(Scenario = rep(1:nrow(settings[[ii]]), each = nrow(uniq_sc[[ii]])), .before = 1) %>%
      select(-c(index)) %>%
      left_join(settings[[ii]], by = "Scenario") %>%
      relocate(names(settings[[ii]]), .after = Scenario)

    output_tbl[[ii]] <- output_tbl[[ii]] %>%
      mutate(Scenario = (floor((Scenario - 1) / 2) + 1))
  }
})


#### Table S3 ===========================================

display_col <- c("Scenario", "Method", "fpr", "fnr", "cfpr_3sigma", "cfpr_2sigma", "cfnr_3sigma", "cfnr_2sigma")
tableS3 <- list()
tableS3_tmp <- list()
rho_source <- 0.5
for (ll in 1:2) {
  beta_shift_type <- c(1, 3)[ll]

  ### iid error
  list_id <- which(sapply(sensitivity_env$output_tbl, function(vv) unique(vv$rho_source)) == rho_source)

  tableS3_tmp[[ll]] <- sensitivity_env$output_tbl[[list_id]] %>%
    filter(beta_delta_idx == beta_shift_type) %>%
    select(all_of(display_col)) %>%
    mutate(
      Method = factor(Method,
        levels = c("onsite_smm", "oracle_smm", "pool_smm", "2stage_smm", "rsmm", "gmm3", "gmm4", "gmm5")
      )
    ) %>%
    mutate(Method = recode(Method,
      "onsite_smm" = "Onsite SMM",
      "oracle_smm" = "Oracle SMM",
      "pool_smm"   = "Pooled SMM",
      "2stage_smm" = "Two-stage SMM",
      "rsmm"       = "Robust SMM",
      "gmm3"       = "Pooled sGMM v0",
      "gmm4"       = "Pooled sGMM",
      "gmm5"       = "TL-sGMM"
    )) %>%
    filter(!(Method %in% c("Onsite SMM", "Oracle SMM", "Pooled sGMM v0")))

  scenario_levels <- sort(unique(tableS3_tmp[[ll]]$Scenario))
  tableS3_tmp[[ll]] <- tableS3_tmp[[ll]] %>%
    mutate(Scenario = factor(Scenario,
      levels = scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
    )) %>%
    mutate(Scenario = recode(
      Scenario,
      !!!setNames(
        c("1", "2a", "2b", "2c", "2d", "3a", "3b"),
        scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
      )
    )) %>%
    arrange(Scenario, Method) %>%
    replace_na(list(fnr = 0, cfnr = 0)) %>%
    mutate(blank1 = NA, .after = 6)

  empty_row <- tableS3_tmp[[ll]][1, ] %>% mutate(across(everything(), ~NA))
  tableS3_tmp[[ll]] <- tableS3_tmp[[ll]] %>%
    group_by(Scenario) %>%
    group_split() %>%
    map_df(~ bind_rows(.x, empty_row)) %>%
    group_by(Scenario) %>%
    mutate(
      Scenario = if_else(row_number() == 1, Scenario, NA_character_)
    ) %>%
    ungroup()
}
tableS3 <- bind_cols(
  tableS3_tmp[[1]],
  NA,
  (tableS3_tmp[[2]] %>% select(-c(Scenario, Method)))
)

#### Table S4 ===========================================
display_col <- c("Scenario", "Method", "rbias_point", "rmse_point", "CI_length", "CI_coverage", "rbias_fwc", "rmse_fwc", "CI_length_fwc", "CI_coverage_fwc")

tableS4 <- list()
tableS4_tmp <- list()
rho_source <- 0.5
for (ll in 1:2) {
  beta_shift_type <- c(1, 3)[ll]
  ### iid error
  list_id <- which(sapply(sensitivity_env$output_tbl, function(vv) unique(vv$rho_source)) == rho_source)
  tableS4_tmp[[ll]] <- sensitivity_env$output_tbl[[list_id]] %>%
    filter(beta_delta_idx == beta_shift_type) %>%
    select(all_of(display_col)) %>%
    mutate(
      Method = factor(Method,
        levels = c("onsite_smm", "oracle_smm", "pool_smm", "2stage_smm", "rsmm", "gmm3", "gmm4", "gmm5")
      )
    ) %>%
    mutate(Method = recode(Method,
      "onsite_smm" = "Onsite SMM",
      "oracle_smm" = "Oracle SMM",
      "pool_smm"   = "Pooled SMM",
      "2stage_smm" = "Two-stage SMM",
      "rsmm"       = "Robust SMM",
      "gmm3"       = "Pooled sGMM v0",
      "gmm4"       = "Pooled sGMM",
      "gmm5"       = "TL-sGMM"
    )) %>%
    filter(!(Method %in% c("Pooled sGMM v0")))

  scenario_levels <- sort(unique(tableS4_tmp[[ll]]$Scenario))
  tableS4_tmp[[ll]] <- tableS4_tmp[[ll]] %>%
    mutate(Scenario = factor(Scenario,
      levels = scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
    )) %>%
    mutate(Scenario = recode(
      Scenario,
      !!!setNames(
        c("1", "2a", "2b", "2c", "2d", "3a", "3b"),
        scenario_levels[c(1, 4, 5, 2, 3, 7, 6)]
      )
    )) %>%
    arrange(Scenario, Method) %>%
    mutate(blank2 = NA, .after = 6)

  empty_row <- tableS4_tmp[[ll]][1, ] %>% mutate(across(everything(), ~NA))
  tableS4_tmp[[ll]] <- tableS4_tmp[[ll]] %>%
    group_by(Scenario) %>%
    group_split() %>%
    map_df(~ bind_rows(.x, empty_row)) %>%
    group_by(Scenario) %>%
    mutate(
      Scenario = if_else(row_number() == 1, Scenario, NA_character_)
    ) %>%
    ungroup()
}

tableS4 <- list(tableS4_tmp[[1]], tableS4_tmp[[2]])



## Paper format: Supplementary Tables S2-S12 ======================================
# percentages to 1 decimal, RMSE and CI length to 2 decimals, halves rounded away from zero
round_half <- function(x, d) {
  r <- sign(x) * floor(abs(x) * 10^d + 0.5 + 1e-8) / 10^d
  r[r == 0] <- 0
  r
}
paper_fmt <- function(df) {
  df <- as.data.frame(df, check.names = FALSE)
  for (j in seq_along(df)) {
    if (!is.numeric(df[[j]]) || names(df)[j] == "L") next
    d <- if (grepl("rmse|CI_length", names(df)[j])) 2 else 1
    df[[j]] <- ifelse(is.na(df[[j]]), NA, sprintf(paste0("%.", d, "f"), round_half(df[[j]], d)))
  }
  df
}
tables_supp <- list(
  "Table S2" = table1[["consistent"]], # tab_detection1
  "Table S3" = table2[["iid error, consistent"]], # tab_est_iid1
  "Table S4" = table2[["iid error, shift"]], # tab_est_iid2
  "Table S5" = table1[["shift"]], # tab_detection2
  "Table S6" = table2[["AR error, consistent"]], # tab_est_ar1
  "Table S7" = table2[["AR error, shift"]], # tab_est_ar2
  "Table S8" = tableS3, # tab_5visits_detect
  "Table S9" = tableS4[[1]], # tab_5visits_est1
  "Table S10" = tableS4[[2]], # tab_5visits_est2
  "Table S11" = tableS1, # tab_sensitivity1
  "Table S12" = tableS2 # tab_sensitivity2
)
write.xlsx(lapply(tables_supp, paper_fmt),
  file = "output/tables/table_simulation.xlsx", asTable = TRUE, withFilter = FALSE
)
