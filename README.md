# Robust Bayesian Transfer Learning and Outlier Detection for Longitudinal Data in Hybrid Decentralized Clinical Trials: An Application to Weight-Loss Trials

This repository holds the code and saved results behind the simulation studies and the real-data application in the paper. Each table in the paper can be regenerated from the saved results in a few minutes, without rerunning any MCMC.

The saved simulation results are too large for this repository and are archived on Zenodo: [doi:10.5281/zenodo.23964967](https://doi.org/10.5281/zenodo.23964967).

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
                          saved results for scenario i (500 replicates);
                          from simulation_results.zip on Zenodo, not in this repository
real_data/
  main.R, stan_more.R     real-data analysis (one model/setting per job; 1,392 settings)
  runjobs.lsf             LSF job array
  make_synthetic_data.R   generates data_synthetic.RData (see "Data availability")
  data_synthetic.RData    synthetic data in the same format as the trial data
  results/output.csv      saved estimates for all 1,392 settings
  extract_shift_posterior.R   extracts shift_posterior.rds from the TL-sGMM fits (not shared)
  results/shift_posterior.rds posterior draws of the shift parameters, used for Figure S1
summary/                  scripts that build the paper's tables and Figure S1
output/                   tables and figures produced by summary/ (as included in this repo)
sessionInfo.txt           R session information for the cluster runs
```

## Requirements

- R 4.3.1. Package versions are listed in `sessionInfo.txt`.
- Stan, through **rstan** 2.32.6.
- JAGS 4.3.x (https://mcmc-jags.sourceforge.io), used by **rjags** and **runjags**.
- The summary scripts also need **openxlsx**, **ggplot2**, **patchwork**, **ggsci**, **ggtext**, **purrr**, **scales** and **RColorBrewer**.

## Reproducing the tables and figures

Run each script from the repository root. Each sheet in the output workbooks is named after its table number in the paper.

**Simulation tables (Tables S2–S12).**

1. Download `simulation_results.zip` from Zenodo ([doi:10.5281/zenodo.23964967](https://doi.org/10.5281/zenodo.23964967)) into the repository root.
2. Unzip it there: `unzip simulation_results.zip`. This places each scenario's result files in its folder under `simulation/`.
3. Run `Rscript summary/tables_simulation.R`. It writes `output/tables/table_simulation.xlsx`:

- Tables S2 and S5: outlier detection with Δ = 0 and Δ = 2
- Tables S3 and S4: estimation with i.i.d. errors, Δ = 0 and Δ = 2
- Tables S6 and S7: estimation with AR errors, Δ = 0 and Δ = 2
- Tables S8–S10: sensitivity analysis, alternative five-visit schedule
- Tables S11 and S12: sensitivity analysis, number of clustering components L

**Application tables (Tables 1, 2 and S1).** Run `Rscript summary/table_application.R`. It writes `output/tables/table_application.xlsx`:

- Table 1: SURMOUNT-OSA treatment effects
- Table 2: TL-sGMM treatment effects with reduced onsite data
- Table S1: arm-specific SURMOUNT-OSA estimates

**Figure S1.** Run `Rscript summary/figure_S1.R` (under a minute). It writes `output/figures/figure_S1.pdf`, the posterior of the TL-sGMM remote–onsite shift for each arm. It reads `real_data/results/shift_posterior.rds`, which holds only the posterior draws of the shift parameters, the spline knots and the arm labels from the four TL-sGMM fits (full onsite data; observed data and induced shift). The full fits contain patient-level data and are not shared; `real_data/extract_shift_posterior.R` shows how the file was extracted from them.

**Figures 1 and 2: not provided.** Both figures show individual patients' data from the SURMOUNT-OSA trial (weight trajectories, and per-patient outlier probabilities), which cannot be shared (see "Data availability"), so they are not reproducible from this repository.

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

The simulation results are archived on Zenodo ([doi:10.5281/zenodo.23964967](https://doi.org/10.5281/zenodo.23964967)). They come from simulated data only.

The real data come from the SURMOUNT-OSA trial and were provided by the industry sponsor under a data use agreement, so they cannot be shared. `real_data/data_synthetic.RData` is fully simulated by `make_synthetic_data.R` and contains no trial information. It has the same objects, column names and types as the trial data, so the real-data analysis code (`real_data/main.R`) can be run end to end. Results obtained from it are not meaningful.
