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

## Tables: Paper ======================================
### Application ==========================
# mirrors code/real_data/main.R: treatment remotes += shift.int + shift.slope * time
# control remotes and all onsite visits are unshifted
# VISIT2 (week 0) dropped from fitting; baseline WEIGHTBL kept as a covariate
# scenarios: (0, 0), (-2, -2); data.type 1:6
# 1 full (post-baseline onsite); 2 weeks 4,12,24,36,52; 3 weeks 4,8,16,24,52
# 4 weeks 4,8,12,36,52; 5 weeks 4,12,36,52; 6 weeks 4,12,52
# MMRM: 12 onsite; 13 pooled (equal-weight onsite + remote within ±5 days of the onsite visit grid)
horseshoe_prior <- data.frame(
  lambda = c(
    "student_t(1, 0, 1)", "student_t(1, 0, 1)", "student_t(1, 0, 1)", "student_t(1, 0, 1)",
    rep("student_t(1, 0, 5)", 2)
  ),
  tau = c(
    "student_t(1, 0, 1)", "student_t(1, 0, 2)", "student_t(1, 0, 5)", "student_t(1, 0, 10)",
    "student_t(1, 0, 5)", "student_t(1, 0, 10)"
  )
)
ig_prior <- data.frame(
  a_slab = c(0.1, 0.001), b_slab = c(0.1, 0.001)
)
all.config <- bind_rows(
  # LMM
  expand.grid(
    trt.group = c(1, 0),
    model.type = c(1:3),
    lmm.type = c(1, 2),
    gmm.type = 1,
    k.knots = c(4),
    data.type = c(1:6),
    sigma_prior = c("student_t(1, 0, 5)"),
    hs_prior = 1,
    ig_prior = 1,
    mcmc_seed = c(123),
    shift.int = c(0, -2),
    shift.slope = c(0, -2)
  ),
  # BLMM
  expand.grid(
    trt.group = c(1, 0),
    model.type = c(4:6),
    lmm.type = c(1, 2),
    gmm.type = 1,
    k.knots = c(4),
    data.type = c(1:6),
    sigma_prior = c("student_t(1, 0, 5)"),
    hs_prior = 1,
    ig_prior = 1,
    mcmc_seed = c(123, 234, 345, 456, 567),
    shift.int = c(0, -2),
    shift.slope = c(0, -2)
  ),
  # sGMM
  expand.grid(
    trt.group = c(1, 0),
    model.type = 7,
    lmm.type = 1,
    gmm.type = c(2, 4, 5),
    k.knots = c(4),
    data.type = c(1:6),
    sigma_prior = c("student_t(1, 0, 5)"),
    hs_prior = 1,
    ig_prior = 1,
    mcmc_seed = c(123, 234, 345, 456, 567),
    shift.int = c(0, -2),
    shift.slope = c(0, -2)
  ),
  # TL-sGMM
  expand.grid(
    trt.group = c(1, 0),
    model.type = 7,
    lmm.type = 1,
    gmm.type = c(11),
    k.knots = c(4),
    data.type = c(1:6),
    sigma_prior = c("student_t(1, 0, 5)"),
    hs_prior = c(1),
    ig_prior = c(2),
    mcmc_seed = c(123, 234, 345, 456, 567),
    shift.int = c(0, -2),
    shift.slope = c(0, -2)
  ),
  # MMRM
  expand.grid(
    trt.group = c(1, 0),
    model.type = c(12:13),
    lmm.type = 1,
    gmm.type = 1,
    k.knots = c(4),
    data.type = c(1:6),
    sigma_prior = c("student_t(1, 0, 5)"),
    hs_prior = 1,
    ig_prior = 1,
    mcmc_seed = c(123),
    shift.int = c(0, -2),
    shift.slope = c(0, -2)
  )
) %>%
  filter((shift.int == 0 & shift.slope == 0) |
           (shift.int == -2 & shift.slope == -2))


data.path <- "real_data/results/"
summary_res <- read_csv(paste0(data.path, "output.csv"), show_col_types = FALSE) %>%
  arrange(V1)
colnames(summary_res) <- c("row_idx", paste0("t = ", 0:12), "auc", "sd_r0", "sd_r1", "converge")
tab_summary <- all.config %>%
  mutate(row_idx = row_number(), .before = 1) %>%
  as_tibble() %>%
  left_join(summary_res, by = "row_idx") %>%
  select(-c(row_idx))
tab_summary <- tab_summary %>%
  mutate(methods = case_when(
    model.type == 1 ~ "Onsite SMM",
    model.type == 2 ~ "Pooled SMM",
    model.type == 3 ~ "Remote SMM",
    model.type == 4 ~ "Onsite BSMM",
    model.type == 5 ~ "Pooled BSMM",
    model.type == 6 ~ "Remote BSMM",
    model.type == 7 & gmm.type == 2 ~ "Onsite sGMM v2",
    model.type == 7 & gmm.type == 4 ~ "Pooled sGMM v2",
    model.type == 7 & gmm.type == 5 ~ "Pooled sGMM v2 heteroRE",
    model.type == 7 & gmm.type == 11 ~ "TL-sGMM",
    model.type == 12 ~ "Onsite MMRM",
    model.type == 13 ~ "Pooled MMRM",
    TRUE ~ NA
  ), .after = 1) %>%
  mutate(data.type = case_when(
    data.type == 1 ~ "Full data",
    data.type == 2 ~ "Reduced data (4, 12, 24, 36, 52)",
    data.type == 3 ~ "Reduced data (4, 8, 16, 24, 52)",
    data.type == 4 ~ "Reduced data (4, 8, 12, 36, 52)",
    data.type == 5 ~ "Reduced data (4, 12, 36, 52)",
    data.type == 6 ~ "Reduced data (4, 12, 52)"
  )) %>%
  select(-c(model.type, gmm.type))

tab_summary <- tab_summary %>%
  mutate(
    tau_prior = horseshoe_prior$tau[hs_prior], .after = "hs_prior"
  ) %>%
  mutate(
    lambda_prior = horseshoe_prior$lambda[hs_prior], .after = "hs_prior"
  )

#### Helpers ==========================
clean_string <- function(x) {
  x %>% str_replace_all("-?\\d+\\.\\d+", ~ as.character(round(as.numeric(.x), 2)))
}
avg_ci <- function(ci_vec) {
  bounds <- str_match(ci_vec, "\\(([^,]+),\\s*([^)]+)\\)")
  lower <- mean(as.numeric(bounds[, 2]), na.rm = TRUE)
  upper <- mean(as.numeric(bounds[, 3]), na.rm = TRUE)
  str_c("(", round(lower, 2), ", ", round(upper, 2), ")")
}
te_from_arms <- function(trt_df, ctrl_df) {
  if (nrow(trt_df) == 0 || nrow(ctrl_df) == 0) {
    return(trt_df[0, , drop = FALSE])
  }
  trt_df %>%
    inner_join(ctrl_df, by = "methods", suffix = c("_t", "_c")) %>%
    mutate(
      `t = 12_est` = `t = 12_est_t` - `t = 12_est_c`,
      `t = 12_se` = sqrt(`t = 12_se_t`^2 + `t = 12_se_c`^2),
      `t = 12_CI` = str_c(
        "(", round(`t = 12_est` - 1.96 * `t = 12_se`, 2), ", ",
        round(`t = 12_est` + 1.96 * `t = 12_se`, 2), ")"
      ),
      auc_est = auc_est_t - auc_est_c,
      auc_se = sqrt(auc_se_t^2 + auc_se_c^2),
      auc_CI = str_c(
        "(", round(auc_est - 1.96 * auc_se, 2), ", ",
        round(auc_est + 1.96 * auc_se, 2), ")"
      ),
      shift.int = shift.int_t,
      shift.slope = shift.slope_t,
      sd_r0 = NA_real_,
      sd_r1 = NA_real_,
      converge = if_else(converge_t == "Yes" & converge_c == "Yes", "Yes", "No"),
      sep = NA
    ) %>%
    select(all_of(names(trt_df)))
}
label_group <- function(df, lab) {
  n <- nrow(df)
  df %>% mutate(
    Group = if (n == 0) NA_character_ else c(lab, rep(NA_character_, n - 1)),
    .before = 1
  )
}

display_col <- c(
  "Group", "Method", "Intercept", "Slope", "sd_r0", "sd_r1", "Converge",
  "Week 52", "Week 52 SE", "Week 52 CI", "sep1", "AUC", "AUC SE", "AUC CI"
)

seed <- c(123, 234, 345, 456, 567)
visit_scen <- list(
  `5v 4-12-24-36-52` = "Reduced data (4, 12, 24, 36, 52)",
  `5v 4-8-16-24-52` = "Reduced data (4, 8, 16, 24, 52)",
  `5v 4-8-12-36-52` = "Reduced data (4, 8, 12, 36, 52)",
  `4v 4-12-36-52` = "Reduced data (4, 12, 36, 52)",
  `3v 4-12-52` = "Reduced data (4, 12, 52)"
)
shift_scen <- list(
  `int0_slope0` = c(int = 0, slope = 0),
  `int-2_slope-2` = c(int = -2, slope = -2)
)

#### Table: tab_app_est (full onsite data) ==========================
table0_full <- list()
shift_tabs <- list()
for (ss_name in names(shift_scen)) {
  ss <- shift_scen[[ss_name]]
  table0_seed <- list()
  for (kk in 1:5) {
    table0_tmp <- list()
    for (ll in 1:2) {
      table0_tmp[[ll]] <- tab_summary %>%
        filter(shift.int == ss[["int"]], shift.slope == ss[["slope"]]) %>%
        filter(!(methods %in% c("Onsite BSMM", "Onsite BSMM (Full)", "Pooled BSMM"))) %>%
        mutate(converge = if_else(is.na(converge) | converge == "", "Yes", converge)) %>%
        filter(trt.group == (ll - 1) & lmm.type %in% c(1, 2) & hs_prior %in% c(1) &
          sigma_prior == "student_t(1, 0, 5)" & mcmc_seed %in% c(123, seed[kk])) %>%
        filter(!(lmm.type == 2 & grepl("SMM", methods))) %>%
        filter(kk == 1 | !(grepl("BSMM", methods) & mcmc_seed == 123)) %>%
        filter(kk == 1 | !(grepl("sGMM", methods) & mcmc_seed == 123)) %>%
        filter(!grepl("Remote", methods, ignore.case = TRUE)) %>%
        filter(data.type == "Full data") %>%
        select(c(methods, shift.int, shift.slope, data.type, `t = 12`, auc, sd_r0, sd_r1, converge)) %>%
        mutate(shift.int = shift.int * (ll == 2), shift.slope = shift.slope * (ll == 2)) %>%
        mutate(methods = factor(methods, levels = c(
          "Onsite SMM", "Pooled SMM",
          "Onsite MMRM", "Pooled MMRM",
          "Onsite sGMM v2",
          "Pooled sGMM v2", "Pooled sGMM v2 heteroRE",
          "TL-sGMM"
        ))) %>%
        arrange(methods) %>%
        mutate(methods = as.character(methods)) %>%
        select(-c(data.type)) %>%
        mutate(
          across(
            c(`t = 12`, auc),
            list(
              est = \(x) as.numeric(str_extract(x, "^-?[\\d.]+")),
              se  = \(x) as.numeric(str_extract(x, "(?<=\\()[^,]+")),
              CI  = \(x) str_c("(", str_extract(x, "(?<=CI: )[^)]+"), ")")
            ),
            .names = "{.col}_{.fn}"
          )
        ) %>%
        select(-`t = 12`, -auc) %>%
        mutate(sep = NA, .before = auc_est)
    }

    table0_te <- te_from_arms(table0_tmp[[2]], table0_tmp[[1]])
    table0_seed[[kk]] <- bind_rows(
      label_group(table0_tmp[[2]], "Treatment"),
      table0_tmp[[2]][NA_integer_, ],
      label_group(table0_tmp[[1]], "Control"),
      table0_tmp[[1]][NA_integer_, ],
      label_group(table0_te, "Treatment effect")
    )
    colnames(table0_seed[[kk]]) <- display_col
  }

  names(table0_seed) <- c("seed1", "seed2", "seed3", "seed4", "seed5")

  shift_tabs[[ss_name]] <- bind_rows(
    lapply(table0_seed, \(df) mutate(df, .row = row_number())),
    .id = "seed"
  ) %>%
    group_by(.row) %>%
    summarise(
      Group = first(na.omit(Group)),
      Method = first(na.omit(Method)),
      Intercept = first(na.omit(Intercept)),
      Slope = first(na.omit(Slope)),
      sd_r0 = median(sd_r0, na.rm = TRUE),
      sd_r1 = median(sd_r1, na.rm = TRUE),
      Converge = paste0(sum(Converge == "Yes", na.rm = TRUE), "/", sum(!is.na(Converge))),
      `Week 52` = mean(`Week 52`, na.rm = TRUE),
      `Week 52 SE` = mean(`Week 52 SE`, na.rm = TRUE),
      `Week 52 CI` = avg_ci(`Week 52 CI`),
      sep1 = first(sep1),
      AUC = mean(AUC, na.rm = TRUE),
      `AUC SE` = mean(`AUC SE`, na.rm = TRUE),
      `AUC CI` = avg_ci(`AUC CI`),
      .groups = "drop"
    ) %>%
    select(-.row)
}
table0_full[["Full data"]] <- bind_rows(
  shift_tabs[[names(shift_scen)[1]]],
  shift_tabs[[names(shift_scen)[1]]][NA_integer_, ],
  shift_tabs[[names(shift_scen)[2]]]
)

#### Table: tab_app_reduced (TL-sGMM, reduced onsite data) ==========================
table0_reduced <- list()
shift_tabs <- list()
for (ss_name in names(shift_scen)) {
  ss <- shift_scen[[ss_name]]
  table0_seed <- list()
  for (kk in 1:5) {
    table0_tmp <- list()
    for (ll in 1:2) {
      table0_tmp[[ll]] <- tab_summary %>%
        filter(shift.int == ss[["int"]], shift.slope == ss[["slope"]]) %>%
        mutate(converge = if_else(is.na(converge) | converge == "", "Yes", converge)) %>%
        filter(trt.group == (ll - 1) & lmm.type %in% c(1, 2) & hs_prior %in% c(1) &
          sigma_prior == "student_t(1, 0, 5)" & mcmc_seed %in% c(123, seed[kk])) %>%
        filter(kk == 1 | !(grepl("sGMM", methods) & mcmc_seed == 123)) %>%
        filter(methods == "TL-sGMM", data.type %in% unlist(visit_scen[c(1, 4, 5)])) %>%
        select(c(methods, shift.int, shift.slope, data.type, `t = 12`, auc, sd_r0, sd_r1, converge)) %>%
        mutate(shift.int = shift.int * (ll == 2), shift.slope = shift.slope * (ll == 2)) %>%
        mutate(methods = factor(data.type, levels = unlist(visit_scen[c(1, 4, 5)]))) %>%
        arrange(methods) %>%
        mutate(methods = as.character(methods)) %>%
        select(-c(data.type)) %>%
        mutate(
          across(
            c(`t = 12`, auc),
            list(
              est = \(x) as.numeric(str_extract(x, "^-?[\\d.]+")),
              se  = \(x) as.numeric(str_extract(x, "(?<=\\()[^,]+")),
              CI  = \(x) str_c("(", str_extract(x, "(?<=CI: )[^)]+"), ")")
            ),
            .names = "{.col}_{.fn}"
          )
        ) %>%
        select(-`t = 12`, -auc) %>%
        mutate(sep = NA, .before = auc_est)
    }

    table0_te <- te_from_arms(table0_tmp[[2]], table0_tmp[[1]])
    table0_seed[[kk]] <- bind_rows(
      label_group(table0_tmp[[2]], "Treatment"),
      table0_tmp[[2]][NA_integer_, ],
      label_group(table0_tmp[[1]], "Control"),
      table0_tmp[[1]][NA_integer_, ],
      label_group(table0_te, "Treatment effect")
    )
    colnames(table0_seed[[kk]]) <- display_col
  }

  names(table0_seed) <- c("seed1", "seed2", "seed3", "seed4", "seed5")

  shift_tabs[[ss_name]] <- bind_rows(
    lapply(table0_seed, \(df) mutate(df, .row = row_number())),
    .id = "seed"
  ) %>%
    group_by(.row) %>%
    summarise(
      Group = first(na.omit(Group)),
      Method = first(na.omit(Method)),
      Intercept = first(na.omit(Intercept)),
      Slope = first(na.omit(Slope)),
      sd_r0 = median(sd_r0, na.rm = TRUE),
      sd_r1 = median(sd_r1, na.rm = TRUE),
      Converge = paste0(sum(Converge == "Yes", na.rm = TRUE), "/", sum(!is.na(Converge))),
      `Week 52` = mean(`Week 52`, na.rm = TRUE),
      `Week 52 SE` = mean(`Week 52 SE`, na.rm = TRUE),
      `Week 52 CI` = avg_ci(`Week 52 CI`),
      sep1 = first(sep1),
      AUC = mean(AUC, na.rm = TRUE),
      `AUC SE` = mean(`AUC SE`, na.rm = TRUE),
      `AUC CI` = avg_ci(`AUC CI`),
      .groups = "drop"
    ) %>%
    select(-.row)
}
table0_reduced[["TL-sGMM"]] <- bind_rows(
  shift_tabs[[names(shift_scen)[1]]],
  shift_tabs[[names(shift_scen)[1]]][NA_integer_, ],
  shift_tabs[[names(shift_scen)[2]]]
)


#### Paper format ==========================
# rows used in the paper, rounded as in the paper
fmt_ci <- function(ci) {
  b <- str_match(ci, "\\(([^,]+),\\s*([^)]+)\\)")
  sprintf("(%.2f, %.2f)", as.numeric(b[, 2]), as.numeric(b[, 3]))
}
paper_rows <- function(df, methods, groups = c("Treatment effect" = "Treatment effect")) {
  df %>%
    fill(Group) %>%
    # each shift scenario is one block starting with the "Treatment" rows
    mutate(Scenario = c("Observed Data", "Induced Shift")[
      cumsum(Group == "Treatment" & lag(Group, default = "") != "Treatment")
    ]) %>%
    filter(Group %in% names(groups), Method %in% names(methods)) %>%
    transmute(
      Scenario,
      Arm = factor(groups[Group], levels = unname(groups)),
      Method = factor(methods[Method], levels = unname(methods)),
      `Week 52 Estimate` = sprintf("%.2f", `Week 52`),
      `Week 52 SE` = sprintf("%.2f", `Week 52 SE`),
      `Week 52 95% CI` = fmt_ci(`Week 52 CI`),
      `AUC Estimate` = sprintf("%.2f", AUC),
      `AUC SE` = sprintf("%.2f", `AUC SE`),
      `AUC 95% CI` = fmt_ci(`AUC CI`)
    ) %>%
    arrange(factor(Scenario, levels = c("Observed Data", "Induced Shift")), Arm, Method) %>%
    mutate(Arm = as.character(Arm), Method = as.character(Method))
}
methods_full <- c(
  "Onsite SMM" = "Onsite SMM", "Pooled SMM" = "Pooled SMM",
  "Onsite MMRM" = "Onsite MMRM", "Pooled MMRM" = "Pooled MMRM", "TL-sGMM" = "TL-sGMM"
)
# Table 1 (tab_app_est): treatment effects, full onsite data
tab_app_est <- paper_rows(table0_full[["Full data"]], methods_full) %>%
  select(-Arm)
# Table 2 (tab_app_reduced): TL-sGMM treatment effects, reduced onsite data
tab_app_reduced <- paper_rows(table0_reduced[["TL-sGMM"]], c(
  "Reduced data (4, 12, 24, 36, 52)" = "5", "Reduced data (4, 12, 36, 52)" = "4",
  "Reduced data (4, 12, 52)" = "3"
)) %>%
  select(-Arm) %>%
  rename(`# Onsite Visits` = Method)
# Table S1 (tab_app_arm): arm-specific estimates, full onsite data
tab_app_arm <- paper_rows(table0_full[["Full data"]], methods_full,
  groups = c("Treatment" = "Tirzepatide", "Control" = "Placebo")
)

write.xlsx(list("Table 1" = tab_app_est, "Table 2" = tab_app_reduced, "Table S1" = tab_app_arm),
  file = "output/tables/table_application.xlsx", asTable = TRUE, withFilter = FALSE
)
