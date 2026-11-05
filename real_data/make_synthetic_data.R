# Generate data_synthetic.RData: a fully simulated dataset of 60 patients (30 per arm)
# with the same objects, column names and types as the (proprietary) trial data
# used in the paper. No values are taken from the trial data; it only illustrates
# the data format expected by main.R and lets every model in it run.
# Results on it are not meaningful.
# Run from the real_data/ folder: Rscript make_synthetic_data.R
library(dplyr)
library(hms)

set.seed(2026)

# patients 1-4: consistent onsite and remote measurements;
# patients 5-8: one remote outlier pattern each; patients 9-60: random patterns
fig_consistent <- 1:4
fig_outlier <- 5:8
n_onsite <- 60
ids_all <- 1:n_onsite
ids_remote <- ids_all

## Baseline characteristics ------------------------------------------------
data_baseline_onsite <- tibble(
  USUBJID = as.numeric(ids_all),
  SSTUDYID = sample(c("GPI1", "GPI2"), n_onsite, replace = TRUE),
  TRT01P = rep(c("TZP MTD", "Placebo"), length.out = n_onsite),
  AGE = as.numeric(sample(30:75, n_onsite, replace = TRUE)),
  SEX = as.numeric(sample(1:2, n_onsite, replace = TRUE)),
  RACE = as.numeric(sample(1:6, n_onsite, replace = TRUE)),
  ETHNIC = as.numeric(sample(1:3, n_onsite, replace = TRUE)),
  ACOUNTRY = as.numeric(sample(1:9, n_onsite, replace = TRUE)),
  HEIGHTBL_perturbed = round(rnorm(n_onsite, 172, 9), 1),
  WEIGHTBL = round(rnorm(n_onsite, 115, 18), 1)
) %>%
  mutate(
    VSSTRESN_perturbed = WEIGHTBL,
    BMIBL_perturbed = round(WEIGHTBL / (HEIGHTBL_perturbed / 100)^2, 1),
    .before = WEIGHTBL
  )
data_baseline <- data_baseline_onsite %>% filter(USUBJID %in% ids_remote)

## Subject-level mean trajectory (weight change, kg) ------------------------
subj <- data_baseline_onsite %>%
  transmute(USUBJID, WEIGHTBL, TRT01P,
    # long-run % weight change: arm mean plus patient-level variation
    target = (if_else(TRT01P == "TZP MTD", -0.20, -0.025) + rnorm(n_onsite, 0, 0.05)) * WEIGHTBL,
    slope_mult = pmax(rnorm(n_onsite, 1, 0.3), 0.2),
    b0 = rnorm(n_onsite, 0, 1.5)
  )
mean_traj <- function(s, week) {
  s$b0 * (week > 0) + s$slope_mult * s$target * (1 - exp(-week / 14))
}

## Onsite visits -------------------------------------------------------------
visits <- tibble(
  VISITNUM = as.numeric(2:11),
  week = c(0, 4, 8, 12, 16, 20, 24, 36, 48, 52)
)
onsite <- subj %>%
  tidyr::crossing(visits) %>%
  filter(VISITNUM == 2 | runif(n()) > 0.08) %>%
  mutate(
    ady = if_else(week == 0, 1, week * 7 + 1 + sample(-3:3, n(), replace = TRUE)),
    time = round((ady - 1) / 30, 3)
  ) %>%
  rowwise() %>%
  mutate(diffw = if_else(week == 0, 0, round(mean_traj(pick(everything()), (ady - 1) / 7) + rnorm(1, 0, 1), 1))) %>%
  ungroup() %>%
  left_join(select(data_baseline_onsite, USUBJID, SSTUDYID), by = "USUBJID") %>%
  transmute(
    SSTUDYID, USUBJID,
    VSGRPID = paste0("VS1001-", VISITNUM, "-1"),
    VISITNUM, VISIT = paste0("VISIT", VISITNUM), ady,
    avisit = paste0("Week ", week, " (Visit ", VISITNUM, ")"),
    VSSTRESN_perturbed = round(WEIGHTBL + diffw, 1),
    hour = NA, time, source = "Onsite", WEIGHTBL, TRT01P, diffw
  )
data_onsite_only <- onsite

## Remote (home) measurements -----------------------------------------------
remote <- lapply(ids_remote, function(id) {
  s <- filter(subj, USUBJID == id)
  n_obs <- min(max(round(exp(rnorm(1, 3, 1.2))), 1), 450)
  if (id %in% c(fig_consistent, fig_outlier)) n_obs <- max(n_obs, 40)
  ady <- sort(sample(2:365, n_obs, replace = n_obs > 364))
  hr <- runif(n_obs, 5, 23)
  week <- (ady - 1 + hr / 24) / 7
  y <- mean_traj(s, week) + rnorm(1, 0, 0.5) + rnorm(n_obs, 0, 1.2)
  # remote-specific outliers
  pattern <- if (id %in% fig_consistent) {
    0
  } else if (id %in% fig_outlier) {
    match(id, fig_outlier)
  } else {
    sample(0:4, 1, prob = c(0.7, 0.1, 0.1, 0.05, 0.05))
  }
  if (pattern == 1) { # sporadic spikes
    k <- sample(n_obs, max(1, round(0.08 * n_obs)))
    y[k] <- y[k] + sample(c(-1, 1), length(k), TRUE) * runif(length(k), 8, 20)
  } else if (pattern == 2) { # persistent shift over a block of time
    k <- which(week > runif(1, 8, 20) & week < runif(1, 30, 45))
    y[k] <- y[k] + runif(1, 5, 10)
  } else if (pattern == 3) { # drift away from the onsite trajectory
    y <- y + 0.25 * week
  } else if (pattern == 4) { # mixture of spikes and a shift
    k <- sample(n_obs, max(1, round(0.05 * n_obs)))
    y[k] <- y[k] - runif(length(k), 8, 20)
    k2 <- which(week > 30)
    y[k2] <- y[k2] + 6
  }
  vnum <- visits$VISITNUM[findInterval(ady, visits$week * 7 + 1)]
  tibble(
    SSTUDYID = data_baseline$SSTUDYID[data_baseline$USUBJID == id],
    USUBJID = id,
    VSGRPID = paste0("REMOTE-", sample(10:12, 1), "-", seq_len(n_obs)),
    VISITNUM = as.numeric(vnum),
    VISIT = if_else(runif(n_obs) < 0.02, NA_character_, paste0("VISIT", vnum)),
    ady = as.numeric(ady),
    hour = hms::hms(hours = floor(hr), minutes = floor((hr %% 1) * 60)),
    time = (ady - 1 + hr / 24) / 30,
    VSSTRESN_perturbed = round(s$WEIGHTBL + y, 1),
    source = "Remote",
    avisit = NA_character_,
    WEIGHTBL = s$WEIGHTBL,
    TRT01P = s$TRT01P,
    diffw = round(y, 1)
  )
}) %>% bind_rows()

data_all <- bind_rows(
  onsite %>%
    filter(USUBJID %in% ids_remote) %>%
    mutate(hour = hms::as_hms(NA_real_)),
  remote
) %>%
  select(
    SSTUDYID, USUBJID, VSGRPID, VISITNUM, VISIT, ady, hour, time,
    VSSTRESN_perturbed, source, avisit, WEIGHTBL, TRT01P, diffw
  ) %>%
  arrange(USUBJID, time)

# the analysis scripts rely on dplyr::select (MASS, loaded indirectly, masks it)
select <- dplyr::select

save(data_all, data_baseline_onsite, data_onsite_only, data_baseline, select,
  file = "data_synthetic.RData"
)
