# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

evolveBO is an R package for Bayesian optimization of Bayesian adaptive clinical trial designs. It implements constrained optimization using heteroskedastic Gaussian process surrogates, constraint-aware acquisition functions (Expected Constrained Improvement), and multi-fidelity simulation budgeting.

## What's New in v0.3.0

Major performance improvements across three implementation phases:

**Phase 1 - Acquisition & Batch Diversity**:
- Improved infeasible region handling via expected constraint violations
- Batch diversity using local penalization (González et al., 2016)
- Numerical stability fixes with epsilon guards
- Impact: 10-20% fewer evaluations to convergence

**Phase 2 - Multi-Fidelity Strategy**:
- Adaptive fidelity selection with cost-aware optimization (new default)
- Three methods: adaptive, staged, threshold
- Custom fidelity costs for non-linear scaling
- Impact: 15-25% better budget utilization

**Phase 3 - Performance Optimizations**:
- Warm-start for GP hyperparameters (30-50% faster surrogate fitting)
- Adaptive candidate pool sizing (scales with dimension)
- Early stopping criterion (saves 10-30% of budget)
- Impact: 30-60% reduction in wall-clock time

**Combined**: 50-70% overall efficiency improvement (expected)

See implementation summaries:
- [PHASE1_IMPLEMENTATION_SUMMARY.md](PHASE1_IMPLEMENTATION_SUMMARY.md)
- [PHASE2_IMPLEMENTATION_SUMMARY.md](PHASE2_IMPLEMENTATION_SUMMARY.md)
- [PHASE3_IMPLEMENTATION_SUMMARY.md](PHASE3_IMPLEMENTATION_SUMMARY.md)
- [ALL_PHASES_COMPLETE.md](ALL_PHASES_COMPLETE.md)

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
   - **Fit surrogates** with warm-start (v0.3.0): `fit_surrogates()` reuses hyperparameters from previous iteration
   - **Generate candidate pool**: Adaptive pool size scales with dimension (500 × d, v0.3.0)
   - **Score candidates**: Acquisition function with improved infeasible handling (v0.3.0)
   - **Select diverse batch** (v0.3.0): Local penalization ensures spatial diversity when `q > 1`
   - **Choose fidelity**: Adaptive cost-aware selection (v0.3.0) or staged/threshold methods
   - **Evaluate and record**: Execute simulations and update history
   - **Check convergence** (v0.3.0): Early stopping if improvement < 0.1% for 5 iterations (configurable via `early_stop` parameter)

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
- **Warm-Start** (v0.3.0): `extract_gp_hyperparams()` extracts lengthscale parameters from previous iteration for faster optimization
- **Prediction**: `predict_surrogates()` returns mean and standard deviation for each metric

### Acquisition Functions (R/acquisition.R)

- **Expected Constrained Improvement (ECI)**: `acq_eci()` computes EI × P(feasible)
- **Feasibility Probability**: Product of univariate Gaussian CDFs for each constraint
- **Improved Infeasible Handling** (v0.3.0): `compute_expected_violation()` guides search toward feasibility boundary using probabilistic constraint violations
- **Numerical Stability** (v0.3.0): Epsilon guards (1e-10) prevent division by zero
- **Batch Diversity** (v0.3.0): `select_batch_local_penalization()` ensures spatially diverse batch selection via local penalization (González et al., 2016)
  - `estimate_lipschitz()`: Extract smoothness constant from GP lengthscales
  - `compute_distances()`: Euclidean distance computation for penalization

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

### Fidelity Selection Methods (v0.3.0)

Three methods available via `fidelity_method` parameter:

**1. Adaptive (Recommended Default)**:
- Cost-aware value-per-cost optimization
- Considers: uncertainty (CV), boundary proximity, acquisition value, budget depletion
- Staged weighting: early iterations favor exploration, late iterations favor exploitation
- Cost exponent increases with iteration and budget usage
- Function: `select_fidelity_adaptive()`

**2. Staged (Simple)**:
- Fixed iteration-based thresholds
- Iterations < 30: low fidelity
- Iterations 30-100: medium fidelity
- Iterations > 100: high fidelity
- Function: `select_fidelity_staged()`

**3. Threshold (Legacy)**:
- Simple feasibility probability thresholds
- P(feasible) ≥ 0.75 → high fidelity
- P(feasible) ≥ 0.40 → medium fidelity
- Otherwise → low fidelity
- Function: `select_fidelity_threshold()`

**Dispatcher**: `select_fidelity_method()` routes to appropriate method

**Custom Costs** (v0.3.0): Use `fidelity_costs` parameter to specify non-linear cost relationships when computational cost doesn't scale linearly with replications (e.g., I/O overhead, parallelization effects)

## Early Stopping Configuration

The `early_stop` parameter controls convergence-based early termination:

```r
bo_calibrate(...,
  early_stop = list(
    enabled = TRUE,      # Enable/disable early stopping
    patience = 5,        # BO iterations without improvement before checking
    threshold = 1e-3,    # Minimum relative improvement (0.1%)
    consecutive = 2      # Consecutive patience windows required to stop
  )
)
```

**Defaults** (if `early_stop = NULL`):
- `enabled = TRUE`
- `patience = 5` iterations
- `threshold = 1e-3` (0.1% relative improvement)
- `consecutive = 2` checks

**To disable**: `early_stop = list(enabled = FALSE)`

**How it works**:
1. After each BO iteration, computes best feasible objective over last `patience` iterations
2. Compares to best objective from before that window
3. If relative improvement < `threshold`, increments counter
4. If counter reaches `consecutive`, optimization stops early
5. Also stops if max acquisition value falls below 1e-6 (all candidates exhausted)

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
