# evolveBO Development Notes

These notes summarise the functionality expected by the manuscript analysis scripts in `/home/naimrashid/Downloads/adaptive-trial-bo-paper/scripts/` so that implementation work can be prioritised accordingly.

## High-priority functions to implement/extend

1. **`bo_calibrate()` (R/bo_calibrate.R)**
   - Replace the placeholder with the full Bayesian optimisation loop.
   - Support heteroskedastic GP surrogates, constraint-aware acquisition (ECI/qEHVI), and multi-fidelity promotion (`R_low`, `R_med`, `R_high`).
   - Return a list containing: `history` (tibble of evaluations), `best_theta`, `surrogates`, `policies`, and `diagnostics` (posterior draws required by sensitivity helpers).

2. **Benchmark helpers (new file, e.g., R/benchmark.R)**
   - `benchmark_methods(sim_fun, bounds, objective, constraints, strategies, bo_args, grid_args, random_args, heuristic_args, simulators_per_eval)` runs repeated calibrations across strategies and returns an S3 object with raw per-run results (strategy, run_id, feasibility, OC estimates, iteration trajectories, simulation budgets).
   - Provide `summarise_benchmark()`, `plot_benchmark_trajectory()`, and `plot_benchmark_efficiency()` for downstream reporting.

3. **Constraint reliability tools (new file, e.g., R/reliability.R)**
   - `estimate_constraint_reliability(...)` reruns strategies across seeds and validates calibrated designs with large Monte Carlo samples.
   - Provide `plot_constraint_reliability()` alongside the summary tibble requested by the scripts.

4. **Multi-fidelity ablation (new file, e.g., R/ablation.R)**
   - `ablation_multifidelity(...)` runs BO under different fidelity policies, capturing simulation budgets and OC performance.
   - `plot_multifidelity_tradeoff()` visualises budget vs accuracy.

5. **Sensitivity diagnostics (extend R/sensitivity.R)**
   - `sensitivity_diagnostics(...)` computes Sobol indices, gradient draws, and kernel/acquisition/prior comparisons.
   - Provide plotting helpers: `plot_sobol_indices()`, `plot_gradient_heatmap()`, `plot_kernel_comparison()`, `plot_acquisition_comparison()`.

6. **Case-study summaries and diagnostics (new file, e.g., R/case_study.R)**
   - `summarise_case_study(fit)` extracts optimal parameters, OC estimates, and credible intervals.
   - `case_study_diagnostics(fit, sobol_samples, gradient_points)` returns Sobol indices, gradient summaries, and gradient covariance matrices for the calibrated design.
   - Plotting helpers: `plot_feasible_frontier(fit)`, `plot_tradeoff_surfaces(fit)`, `plot_case_sobol()`, `plot_case_gradient()`, `plot_case_covariance()` so the manuscript scripts can render the Section 4 figures.

## Package infrastructure

- Accept a `sim_fun` callback with signature `function(theta, fidelity = c("low","med","high"), ...)` (matches `sim_fun_evolveTrial()` in the manuscript repo).
- Add roxygen documentation and export tags for new functions; regenerate `NAMESPACE` with roxygen2.
- Add smoke tests under `tests/testthat/` covering benchmarking, reliability, and sensitivity helpers.
- Consider vignettes demonstrating end-to-end usage with the adaptive-trial example.

## Data structures expected by scripts

- Benchmark object should store:
  ```
  list(
    results = tibble(strategy, run_id, feasible, oc_power, oc_type1, oc_en, oc_et,
                     total_sim_calls, history = list(tibble(iter, objective, feasible, fidelity))),
    summary = tibble(...)
  )
  ```
- Reliability object: `summary` (strategy-level feasibility rates) and `runs` tibble with per-run outcomes.
- Multi-fidelity object: `summary` with columns `policy`, `sim_calls`, `objective`, `constraint_gap`.
- Sensitivity object: list with `sobol`, `gradients`, `kernel_comparison`, `acquisition_comparison`.
- Case-study summary: `design` tibble (parameter, value) and `operating_characteristics` tibble (metric, estimate, lower, upper).

## External dependencies

- Likely packages: `hetGP`, `DiceKriging`, `lhs`, `tgp`/`sobolx`, `tidyverse`, `furrr`/`future.apply` for parallel loops.
- Update DESCRIPTION Imports/Suggests as needed.

## Open questions

- Define a reproducible heuristic tuning procedure (e.g., coarse coordinate search).
- Confirm simulator fidelity levels align with Section 3 defaults (200/1000/10000 replicates).
- Ensure plotting helpers produce publication-ready figures (consistent theme, labelled axes, legend text).

Once these functions are implemented, the analysis scripts in the manuscript repository will run end-to-end.
