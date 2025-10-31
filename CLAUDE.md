# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

evolveBO is an R package for Bayesian optimization of Bayesian adaptive clinical trial designs. It implements constrained optimization using heteroskedastic Gaussian process surrogates, constraint-aware acquisition functions (Expected Constrained Improvement), and multi-fidelity simulation budgeting.

## Development Commands

### Building and Documentation
```r
# Generate documentation from roxygen comments
roxygen2::roxygenise()

# Build package
R CMD build .

# Install from source
R CMD INSTALL .

# Or install with devtools
devtools::install()
devtools::load_all()  # Load package for development
```

### Testing
```r
# Run all tests
testthat::test_check("evolveBO")

# Run tests interactively
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-evolveBO-core.R")
```

### Checking Package
```r
# Full R CMD check
devtools::check()

# Quick check without examples/vignettes
devtools::check(document = FALSE, vignettes = FALSE)
```

## Core Architecture

### Simulator Contract
All calibration functions accept a `sim_fun` callback with signature:
```r
function(theta, fidelity = c("low", "med", "high"), seed = NULL, ...)
```
- **Input**: `theta` is a named list/vector of design parameters
- **Output**: Named numeric vector of operating characteristics (e.g., power, type1, EN, ET)
- **Attributes**: Optionally attach `variance` (Monte Carlo variances) and `n_rep` (replication count)

### Main Optimization Pipeline (R/bo_calibrate.R)

The `bo_calibrate()` function orchestrates the complete Bayesian optimization loop:

1. **Initial Design**: Latin hypercube sampling (`lhs_design()`) generates `n_init` starting points
2. **History Tracking**: Each evaluation stores `theta`, `unit_x` (scaled to [0,1]), metrics, variance, feasibility
3. **Iterative Loop**:
   - Fit heteroskedastic GP surrogates for objective + constraints (`fit_surrogates()`)
   - Generate candidate pool via LHS (`lhs_candidate_pool()`)
   - Score candidates using acquisition function (`evaluate_acquisition()`)
   - Select top `q` candidates for evaluation
   - Dynamically choose fidelity based on predicted feasibility (`select_fidelity()`)
   - Record new evaluations and repeat

**Returns**: `evolveBO_fit` object containing:
- `history`: tibble of all evaluations
- `best_theta`: optimal design parameters
- `surrogates`: fitted GP models
- `policies`: configuration settings
- `diagnostics`: posterior draws for sensitivity analysis

### Surrogate Modeling (R/surrogates.R)

- **Aggregation**: Evaluations at identical `theta` (via `theta_id`) are averaged before fitting
- **Noise Handling**: If simulator returns variance estimates, uses heteroskedastic GP; otherwise adds small nugget
- **DiceKriging**: Uses `km()` with Matérn kernels (default 5/2)
- **Prediction**: `predict_surrogates()` returns mean and standard deviation for each metric

### Acquisition Functions (R/acquisition.R)

- **Expected Constrained Improvement (ECI)**: `acq_eci()` computes EI × P(feasible)
- **Feasibility Probability**: Product of univariate Gaussian CDFs for each constraint
- **Handling Infeasibility**: If no feasible points yet, acquisition favors exploration (standard deviation)

### Constraints (R/constraints.R)

- **Format**: Named list mapping metric → `c(direction, threshold)`
  - `direction`: "ge" (greater/equal) or "le" (less/equal)
- **Deterministic Check**: `is_feasible()` checks if observed metrics satisfy all constraints
- **Probabilistic Check**: `prob_feasibility()` computes P(feasible) under GP posterior

### Parameter Transformations (R/utils.R)

- **Unit Hypercube**: Internally operates on [0,1]^d via `scale_to_unit()` / `scale_from_unit()`
- **Integer Parameters**: `integer_params` argument rounds specified parameters before simulation
- **Theta Hashing**: `theta_to_id()` creates unique string identifier for duplicate detection

## Supporting Modules

The package includes additional analysis helpers (see NOTES_FOR_CODEX.md for implementation details):

- **R/benchmark.R**: Compare BO vs random/grid/heuristic strategies
- **R/reliability.R**: Validate constraint satisfaction with large-sample simulation
- **R/ablation.R**: Multi-fidelity policy ablation studies
- **R/sensitivity.R**: Sobol indices, gradient analysis, kernel/acquisition comparisons
- **R/case_study.R**: Extract and visualize calibrated designs

## Key Data Structures

### History Tibble
Each row represents one simulator evaluation:
- `iter`: BO iteration (0 = initial design)
- `eval_id`: Global evaluation counter
- `theta`: Named list of raw parameter values
- `unit_x`: Named numeric vector scaled to [0,1]
- `theta_id`: String hash for duplicate detection
- `fidelity`: "low"/"med"/"high"
- `n_rep`: Actual simulator replications used
- `metrics`: Named list of operating characteristics
- `variance`: Named list of Monte Carlo variances
- `objective`: Value of objective metric
- `feasible`: Boolean constraint satisfaction

### Constraint Table
Output of `parse_constraints()`:
- `metric`: Name of operating characteristic
- `direction`: "ge" or "le"
- `threshold`: Numeric bound

## Variance Estimation Utilities

**NEW**: The package now provides helper functions for memory-efficient variance computation:

### `welford_mean_var()`
Computes mean and variance incrementally using Welford's algorithm:
- **Memory**: O(m) where m = number of metrics (typically ~64 bytes)
- **Overhead**: ~2% computational cost
- **Benefit**: Enables heteroskedastic GP (30-50% better BO performance)

**Usage**:
```r
sim_fun <- function(theta, fidelity = "high", ...) {
  result <- welford_mean_var(
    sample_fn = function(i, theta) run_one_trial(theta),
    n_samples = 1000,
    theta = theta
  )

  metrics <- result$mean
  attr(metrics, "variance") <- result$variance
  return(metrics)
}
```

### `pool_welford_results()`
Pools variance estimates from parallel simulation chunks using Chan's algorithm.

**See**: `inst/examples/simulator_with_variance.R` for complete examples.

## Dependencies

Core dependencies (see DESCRIPTION):
- **Surrogate modeling**: DiceKriging, hetGP
- **Design of experiments**: lhs, randtoolbox
- **Optimization**: DiceOptim
- **Data manipulation**: tidyverse (dplyr, purrr, tibble, tidyr)
- **Parallelization**: furrr, progressr
- **Statistics**: mvtnorm

## Fidelity Levels

Default multi-fidelity policy (customizable via `fidelity_levels` argument):
- `low`: 200 simulator replications
- `med`: 1000 replications
- `high`: 10000 replications

Fidelity selection heuristic in `select_fidelity()`:
- P(feasible) ≥ 0.75 → high fidelity
- P(feasible) ≥ 0.40 → medium fidelity
- Otherwise → low fidelity

## Common Development Patterns

### Adding a New Acquisition Function
1. Implement function in R/acquisition.R following signature of `acq_eci()`
2. Add case to `evaluate_acquisition()` in R/bo_calibrate.R
3. Update `acquisition` argument choices in `bo_calibrate()` documentation
4. Export function in NAMESPACE via roxygen `@export` tag

### Adding a New Surrogate Type
1. Extend `fit_surrogates()` to accept new `surrogate_type` argument
2. Ensure `predict_surrogates()` handles new model class
3. Update tests to cover new surrogate option

### Testing Workflow
- Use `toy_sim_fun` in tests/testthat/test-evolveBO-core.R as template for simple test cases
- Keep test budgets small (budget ≤ 10) to ensure fast execution
- Mark tests with `skip_on_cran()` if they require external packages or are slow
