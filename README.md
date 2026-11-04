# Robust Bayesian Transfer Learning and Outlier Detection for Longitudinal Data in Hybrid Decentralized Clinical Trials: An Application to Weight-Loss Trials

This repository holds the code and saved results behind the simulation studies and the real-data application in the paper. Each table in the paper can be regenerated from the saved results in a few minutes, without rerunning any MCMC.

## Repository layout

```
simulation/
  iid_error/              main simulations, i.i.d. errors       (352 scenarios)
  ar_error/               main simulations, AR(1) errors        (336 scenarios)
  sensitivity/
    num_clusters/         sensitivity: number of mixture components L   (42 scenarios)
    num_measures/         sensitivity: alternative five-visit schedule  (224 scenarios)
    each folder:
      settings.R          scenario grid (all.config); row i = scenario i
      main.R              runs one replicate of one scenario
      functions.R, stan_functions.R   model fitting functions and Stan code
      submit_jobs.sh, runjobs.lsf     LSF job submission
      output_<i>.csv, muhat_fixed_<i>.csv, mutrue_fixed_<i>.csv
                          saved results for scenario i (500 replicates)
real_data/
  main.R, stan_more.R     real-data analysis (one model/setting per job; 1,392 settings)
  runjobs.lsf             LSF job array
  make_synthetic_data.R   generates data_synthetic.RData (see "Data availability")
  data_synthetic.RData    synthetic data in the same format as the trial data
  results/output.csv      saved estimates for all 1,392 settings
summary/                  scripts that build the paper's tables and figures
output/tables/            tables produced by summary/ (as included in this repo)
sessionInfo.txt           R session information for the cluster runs
```

## Requirements

- R 4.3.1. Package versions are listed in `sessionInfo.txt`.
- Stan, through **rstan** 2.32.6.
- JAGS 4.3.x (https://mcmc-jags.sourceforge.io), used by **rjags** and **runjags**.
- The summary scripts also need **openxlsx**, **ggplot2**, **patchwork**, **ggsci**, **ggtext**, **purrr**, **scales** and **RColorBrewer**.

## Reproducing the tables and figures (minutes)

Run each script from the repository root:

| Paper output | Script | Output file (`output/`) |
|---|---|---|
| Tables S2, S5: outlier detection, Δ = 0 / Δ = 2 (`tab_detection1`, `tab_detection2`) | `Rscript summary/tables_simulation.R` | `tables/table_simulation.xlsx` (sheets `Table S2`, `Table S5`) |
| Tables S3, S4: estimation, i.i.d. errors, Δ = 0 / Δ = 2 (`tab_est_iid1`, `tab_est_iid2`) | same script | `tables/table_simulation.xlsx` (sheets `Table S3`, `Table S4`) |
| Tables S6, S7: estimation, AR errors, Δ = 0 / Δ = 2 (`tab_est_ar1`, `tab_est_ar2`) | same script | `tables/table_simulation.xlsx` (sheets `Table S6`, `Table S7`) |
| Tables S8–S10: alternative five-visit schedule (`tab_5visits_detect`, `tab_5visits_est1`, `tab_5visits_est2`) | same script | `tables/table_simulation.xlsx` (sheets `Table S8`–`Table S10`) |
| Tables S11, S12: number of clustering components L (`tab_sensitivity1`, `tab_sensitivity2`) | same script | `tables/table_simulation.xlsx` (sheets `Table S11`, `Table S12`) |
| Table 1: SURMOUNT-OSA treatment effects (`tab_app_est`) | `Rscript summary/table_application.R` | `tables/table_application.xlsx` (sheet `Table 1`) |
| Table 2: TL-sGMM with reduced onsite data (`tab_app_reduced`) | same script | `tables/table_application.xlsx` (sheet `Table 2`) |
| Table S1: arm-specific SURMOUNT-OSA estimates (`tab_app_arm`) | same script | `tables/table_application.xlsx` (sheet `Table S1`) |
| Figure 1: weight-change trajectories | `Rscript summary/figures_application.R` | `figures/figure1.pdf` * |
| Figure 2: detected outliers, representative subjects | same script | `figures/figure2.pdf` ** |
| Figures 2b, 2c: outlier probabilities for all subjects; posterior of the remote–onsite shift | same script | `figures/figure2b.pdf`, `figures/figure2c.pdf` ** |

`tables_simulation.R` takes about 2–3 minutes. The other scripts take under a minute.

\* With the public synthetic data, Figure 1 shows synthetic trajectories, not the trial data.
\*\* These figures need the posterior fits `real_data/results/fit_{1225,1226,1285,1286}.RData`. These are the TL-sGMM fits with full onsite data (MCMC seed 123), one per arm, without (1225/1226) and with (1285/1286) the induced shift. The fits contain patient-level trial data, so they are not included, and the figures cannot be reproduced from this repository. With access to the trial data, run real-data settings 1225, 1226, 1285 and 1286 (see below) and then `summary/figures_application.R`. Running the same four settings on `data_synthetic.RData` checks that the code runs end to end.

## Rerunning the analyses (HPC, optional)

The saved results above were produced on an LSF cluster. Each simulation replicate and each real-data setting ran as an independent job with 3 cores (one per MCMC chain) and up to 10 GB of memory. Depending on the scenario, a job took roughly 1–30 hours.

**Simulations.** Inside a simulation folder (e.g. `simulation/iid_error/`):

1. Move the saved `output_*.csv`, `muhat_fixed_*.csv` and `mutrue_fixed_*.csv` files out of the folder.
2. In `runjobs.lsf`, replace `/path/to/repo` and `/path/to/logs`, and adjust the `module load` lines for your cluster.
3. Run `bash submit_jobs.sh`. This creates `jobtype<i>/` for each scenario and submits 500 replicates per scenario as LSF job arrays. Each job appends its replicate to `jobtype<i>/results/*_<i>.csv`.
4. When all jobs have finished, move the results up to the folder that the summary scripts read: `mv jobtype*/results/*.csv .`

To run a single replicate without LSF, from inside `jobtype<i>/`:

```bash
perl -pe 's/rep00/1/g; s/seed00/12331/g' main.R > main_1.R && Rscript main_1.R
```

`submit_jobs.sh` has already replaced `sc00` with `<i>` in that folder's `main.R`.

**Real data.** In `real_data/`, `runjobs.lsf` submits settings 1–1,392. To run a single setting `k`:

```bash
perl -pe 's/sc00/k/g' main.R > main_k.R && Rscript main_k.R
```

This writes `results/fit_k.RData` and appends one row to `results/output.csv`. Remove the saved `results/output.csv` before a full rerun.

## Data availability

The real data come from the SURMOUNT-OSA trial and were provided by the industry sponsor under a data use agreement, so they cannot be shared. `real_data/data_synthetic.RData` is fully simulated by `make_synthetic_data.R` and contains no trial information. It has the same objects, column names and types as the trial data, so `real_data/main.R` and the figure scripts run end to end. Results obtained from it are not meaningful.
