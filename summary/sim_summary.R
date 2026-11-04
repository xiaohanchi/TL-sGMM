# Summarise one simulation experiment: reads <result.path>/{output,muhat_fixed,mutrue_fixed}_<i>.csv
# for every row of all.config (from that experiment's settings.R) and builds
# output / output_methods / settings / uniq_sc. Sourced by tables_simulation.R.
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
options(nwarnings = 10000)
rerun_idx <- c()
# Common converged replicates across the 7 paper methods for each DGM scenario
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
  if (n_all < 495) {
    warning(paste0("File output_", ii, " has only ", n_all, " rows (< 500)"))
    rerun_idx <- c(rerun_idx, ii)
  }
  fpr <- round(mean(tmp_res$fpr) * 100, 2)
  fnr <- round(mean(tmp_res$fnr) * 100, 2)
  cfpr_3sigma <- round(mean(tmp_res$cfpr_3sigma) * 100, 2)
  cfnr_3sigma <- round(mean(tmp_res$cfnr_3sigma) * 100, 2)
  cfpr_2sigma <- round(mean(tmp_res$cfpr_2sigma) * 100, 2)
  cfnr_2sigma <- round(mean(tmp_res$cfnr_2sigma) * 100, 2)
  # fpr_pts <- round(mean(tmp_res$fpr_p_num/tmp_res$fpr_p_denom) * 100, 2)
  # if(all(tmp_res$fnr_p_denom > 0)) {
  #   fnr_pts <- round(mean(tmp_res$fnr_p_num/tmp_res$fnr_p_denom) * 100, 2)
  # } else {
  #   fnr_pts <- NA
  # }

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

settings <- lapply(
  seq_along(unique(all.config$rho_source)), function(rr) {
    all.config[uniq_sc[[rr]][1, ], col.idx] %>%
      # group_by(re_factor, beta_delta_idx) %>%
      mutate(Scenario = 1:n(), .before = 1)
  }
  # ungroup()
)
