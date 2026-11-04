# Extract the group-level quantities needed for Figure S1 from the TL-sGMM fits.
# The fits (results/fit_<k>.RData, written by main.R) also hold patient-level data
# (df.fit), so they are not shared; only the output of this script is.
# Run from the real_data/ folder: Rscript extract_shift_posterior.R
library(rstan)

# TL-sGMM, full onsite data, MCMC seed 123 (rows of all.config in main.R)
fits <- c(
  observed_tirzepatide = 1225, observed_placebo = 1226,
  induced_shift_tirzepatide = 1285, induced_shift_placebo = 1286
)

shift_posterior <- lapply(fits, function(k) {
  res <- readRDS(paste0("results/fit_", k, ".RData"))
  db <- rstan::extract(res$stan_fit, pars = "delta_beta")$delta_beta
  list(
    delta_beta = db,
    # interior spline knots, as in set_splines() in main.R
    knots = quantile(unique(res$df.fit$time), seq(0, 1, length = ncol(db))[-c(1, ncol(db))]),
    treatment = res$df.fit$treatment[1]
  )
})

saveRDS(shift_posterior, "results/shift_posterior.rds")
