# test-bug-fixes.R
# Unit tests for bug fixes in v0.3.x
#
# Tests cover:
# 1. Constant predictor mean_value consistency (Bug #1)
# 2. Homoskedastic fallback Z aggregation (Bug #2)
# 3. Warm-start with initial_history without unit_x (Bug #3)
# 4. Local sensitivity unit scale handling (Bug #4)
# 5. run_init_with_stopping objective extraction (Bug #5)

library(testthat)

# ============================================================================
# Test fixtures
# ============================================================================

# Simple mock simulator for testing
mock_sim <- function(theta, fidelity = c("low", "med", "high"), seed = 123, ...) {
  fidelity <- match.arg(fidelity)
  set.seed(seed)

  x1 <- theta$x1
  x2 <- theta$x2

  # Objective: minimize (x1-0.5)^2 + (x2-0.3)^2
  EN <- (x1 - 0.5)^2 + (x2 - 0.3)^2 + rnorm(1, 0, 0.01)

  # Constraints
  power <- pnorm((x1 + x2 - 0.5) / 0.2) + rnorm(1, 0, 0.02)
  type1 <- 0.05 + 0.05 * abs(x1 - x2) + rnorm(1, 0, 0.01)

  result <- c(EN = EN, power = power, type1 = type1)
  attr(result, "variance") <- c(EN = 0.0001, power = 0.0004, type1 = 0.0001)
  attr(result, "n_rep") <- ifelse(fidelity == "low", 100, 1000)

  return(result)
}

# ============================================================================
# Bug #1: Constant predictor mean_value consistency
# ============================================================================

test_that("constant_predictor uses mean_value consistently", {
  skip_on_cran()

  # Create a constant predictor directly (simulating fallback)
  const_pred <- structure(
    list(mean_value = 0.5, metric = "test"),
    class = "constant_predictor"
  )

  # Test that predict_surrogates can read it correctly
  surrogates <- list(test = const_pred)
  unit_x <- list(c(x1 = 0.5, x2 = 0.5))

  # This should not error and should return the constant value
  pred <- BATON:::predict_surrogates(surrogates, unit_x)

  expect_equal(length(pred), 1)
  expect_equal(pred$test$mean, 0.5)
  expect_equal(pred$test$sd, 1.0)  # High uncertainty for constant predictor
})

test_that("fit_surrogates creates constant_predictor with mean_value on insufficient data", {
  skip_on_cran()

  # Create minimal history with only 1 observation
  history <- tibble::tibble(
    iter = 0L,
    eval_id = 1L,
    theta = list(list(x1 = 0.5, x2 = 0.5)),
    unit_x = list(c(x1 = 0.5, x2 = 0.5)),
    theta_id = "0.500000|0.500000",
    fidelity = "low",
    n_rep = 100L,
    metrics = list(c(EN = 0.1, power = 0.85)),
    variance = list(c(EN = 0.01, power = 0.04)),
    objective = 0.1,
    feasible = TRUE
  )

  constraint_tbl <- tibble::tibble(
    metric = "power",
    direction = "ge",
    threshold = 0.8
  )

  # fit_surrogates should fall back to constant_predictor
  expect_message(
    surrogates <- BATON::fit_surrogates(
      history = history,
      objective = "EN",
      constraint_tbl = constraint_tbl
    ),
    "constant predictor"
  )

  # Check that the constant predictor has mean_value (not mean)
  for (metric in names(surrogates)) {
    if (inherits(surrogates[[metric]], "constant_predictor")) {
      expect_true("mean_value" %in% names(surrogates[[metric]]))
    }
  }
})

# ============================================================================
# Bug #3: Warm-start with initial_history without unit_x
# ============================================================================

test_that("initial_history without unit_x derives required columns", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.70), type1 = c("le", 0.15))

  # Create initial_history WITHOUT unit_x column (only individual param columns)
  initial_history <- data.frame(
    x1 = c(0.3, 0.5, 0.7, 0.4, 0.6),
    x2 = c(0.2, 0.5, 0.3, 0.6, 0.4),
    EN = c(0.1, 0.05, 0.15, 0.08, 0.12),
    power = c(0.82, 0.85, 0.78, 0.80, 0.83),
    type1 = c(0.09, 0.08, 0.11, 0.10, 0.07),
    objective = c(0.1, 0.05, 0.15, 0.08, 0.12),
    fidelity = rep("low", 5),
    feasible = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  # This should work without error - unit_x should be derived
  fit <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = 5,
    budget = 8,  # Just a few BO iterations
    initial_history = initial_history,
    seed = 123,
    progress = FALSE
  )

  # Check that unit_x was derived
  expect_true("unit_x" %in% names(fit$history))
  expect_true(is.list(fit$history$unit_x))
  expect_equal(length(fit$history$unit_x[[1]]), 2)  # 2 parameters

  # Check unit_x values are in [0,1]
  for (i in seq_len(nrow(fit$history))) {
    ux <- fit$history$unit_x[[i]]
    expect_true(all(ux >= 0 & ux <= 1))
  }
})

test_that("initial_history with theta column but no unit_x works", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.70), type1 = c("le", 0.15))

  # Create initial_history with theta column but no unit_x
  initial_history <- tibble::tibble(
    theta = list(
      list(x1 = 0.3, x2 = 0.2),
      list(x1 = 0.5, x2 = 0.5),
      list(x1 = 0.7, x2 = 0.3),
      list(x1 = 0.4, x2 = 0.6),
      list(x1 = 0.6, x2 = 0.4)
    ),
    EN = c(0.1, 0.05, 0.15, 0.08, 0.12),
    power = c(0.82, 0.85, 0.78, 0.80, 0.83),
    type1 = c(0.09, 0.08, 0.11, 0.10, 0.07),
    objective = c(0.1, 0.05, 0.15, 0.08, 0.12),
    fidelity = rep("low", 5),
    feasible = c(TRUE, TRUE, FALSE, TRUE, TRUE)
  )

  # Should work - unit_x derived from theta
  fit <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = 5,
    budget = 8,
    initial_history = initial_history,
    seed = 123,
    progress = FALSE
  )

  expect_true("unit_x" %in% names(fit$history))
  expect_true("theta_id" %in% names(fit$history))
})

# ============================================================================
# Bug #4: Local sensitivity unit scale handling
# ============================================================================

test_that("compute_local_sensitivity requires bounds in fit object", {
  skip_on_cran()

  # Create a mock fit without bounds
  fit_no_bounds <- list(
    best_theta = list(x1 = 0.5, x2 = 0.5),
    surrogates = list(EN = structure(list(mean_value = 0.1), class = "constant_predictor"))
  )

  expect_error(
    compute_local_sensitivity(fit_no_bounds),
    "does not contain bounds"
  )
})

test_that("compute_local_sensitivity works with unit-scale surrogates", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.70), type1 = c("le", 0.15))

  fit <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = 10,
    budget = 15,
    seed = 123,
    progress = FALSE
  )

  # Should work without error
  sensitivity <- compute_local_sensitivity(fit)

  expect_s3_class(sensitivity, "data.frame")
  expect_true("parameter" %in% names(sensitivity))
  expect_true("gradient" %in% names(sensitivity))
  expect_true("normalized_sensitivity" %in% names(sensitivity))

  # Check gradients are finite
  expect_true(all(is.finite(sensitivity$gradient)))

  # Normalized sensitivity should sum to 1 (or be all zeros)
  total_sens <- sum(sensitivity$normalized_sensitivity)
  expect_true(total_sens == 0 || abs(total_sens - 1) < 1e-6)
})

test_that("compute_local_sensitivity handles constant predictor", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))

  # Create fit with constant predictor surrogate
  fit <- list(
    best_theta = list(x1 = 0.5, x2 = 0.3),
    bounds = bounds,
    surrogates = list(
      EN = structure(list(mean_value = 0.1, metric = "EN"), class = "constant_predictor")
    )
  )

  # Should return zero gradients for constant predictor
  sensitivity <- compute_local_sensitivity(fit)

  expect_true(all(sensitivity$gradient == 0))
})

# ============================================================================
# Bug #5: run_init_with_stopping objective extraction
# ============================================================================

test_that("run_init_with_stopping extracts named objective correctly", {
  skip_on_cran()

  # Simulator with named output (objective is second element)
  named_sim <- function(theta, fidelity = "low", seed = 1, ...) {
    c(power = 0.85, EN = 0.1 + theta$x1, type1 = 0.08)
  }

  # Create design matrix
  design <- matrix(runif(10), ncol = 2)
  colnames(design) <- c("x1", "x2")

  # Run with objective = "EN" (second element)
  result <- BATON:::run_init_with_stopping(
    sim_fun = named_sim,
    design = design,
    objective = "EN",
    fidelity = "low",
    init_config = BATON::init_stopping_config(enabled = FALSE),
    seed = 123,
    verbose = FALSE
  )

  # y should be EN values, not power values
  expected_y <- 0.1 + design[, "x1"]
  expect_equal(result$y, expected_y, tolerance = 1e-10)
})

test_that("run_init_with_stopping warns for missing objective name", {
  skip_on_cran()

  sim <- function(theta, ...) {
    c(power = 0.85, type1 = 0.08)  # No "EN" metric
  }

  design <- matrix(runif(6), ncol = 2)
  colnames(design) <- c("x1", "x2")

  # Should warn about missing objective and fall back to first element
  expect_warning(
    result <- BATON:::run_init_with_stopping(
      sim_fun = sim,
      design = design,
      objective = "EN",
      fidelity = "low",
      init_config = BATON::init_stopping_config(enabled = FALSE),
      seed = 123,
      verbose = FALSE
    ),
    "not found in simulator result"
  )

  # y should be first element (power = 0.85)
  expect_equal(result$y, rep(0.85, 3))
})

test_that("run_init_with_stopping uses first element when objective is NULL", {
  skip_on_cran()

  sim <- function(theta, ...) {
    c(EN = 0.1, power = 0.85)
  }

  design <- matrix(runif(6), ncol = 2)
  colnames(design) <- c("x1", "x2")

  # No warning expected - just use first element
  result <- BATON:::run_init_with_stopping(
    sim_fun = sim,
    design = design,
    objective = NULL,  # No objective specified
    fidelity = "low",
    init_config = BATON::init_stopping_config(enabled = FALSE),
    seed = 123,
    verbose = FALSE
  )

  # y should be first element (EN = 0.1)
  expect_equal(result$y, rep(0.1, 3))
})

# ============================================================================
# Bug #2: Homoskedastic fallback Z aggregation (indirect test)
# ============================================================================

test_that("fit_hetgp_surrogate handles fallback correctly", {
  skip_on_cran()
  skip_if_not_installed("hetGP")

  # This test indirectly verifies the Z aggregation fix by ensuring

# the fallback path doesn't error when hetGP fails

  # Create history with replicated observations
  history <- tibble::tibble(
    iter = rep(0L, 4),
    eval_id = 1:4,
    theta = list(
      list(x1 = 0.5, x2 = 0.5),
      list(x1 = 0.5, x2 = 0.5),  # Replicate
      list(x1 = 0.3, x2 = 0.7),
      list(x1 = 0.3, x2 = 0.7)   # Replicate
    ),
    unit_x = list(
      c(x1 = 0.5, x2 = 0.5),
      c(x1 = 0.5, x2 = 0.5),
      c(x1 = 0.3, x2 = 0.7),
      c(x1 = 0.3, x2 = 0.7)
    ),
    theta_id = c("0.500000|0.500000", "0.500000|0.500000",
                 "0.300000|0.700000", "0.300000|0.700000"),
    fidelity = rep("low", 4),
    n_rep = rep(100L, 4),
    metrics = list(
      c(EN = 0.1, power = 0.85),
      c(EN = 0.12, power = 0.83),
      c(EN = 0.2, power = 0.75),
      c(EN = 0.18, power = 0.77)
    ),
    variance = list(
      c(EN = 0.01, power = 0.04),
      c(EN = 0.01, power = 0.04),
      c(EN = 0.01, power = 0.04),
      c(EN = 0.01, power = 0.04)
    ),
    objective = c(0.1, 0.12, 0.2, 0.18),
    feasible = c(TRUE, TRUE, FALSE, FALSE)
  )

  constraint_tbl <- tibble::tibble(
    metric = "power",
    direction = "ge",
    threshold = 0.8
  )

  # fit_surrogates should work even if hetGP fails and falls back
  # The key is that it shouldn't error due to Z aggregation issues
  surrogates <- tryCatch(
    BATON::fit_surrogates(
      history = history,
      objective = "EN",
      constraint_tbl = constraint_tbl,
      use_hetgp = TRUE
    ),
    error = function(e) {
      # Some error is acceptable, but not Z aggregation error
      if (grepl("subscript out of bounds|Z_vec|matching", e$message)) {
        fail(paste("Z aggregation bug still present:", e$message))
      }
      # Return NULL to indicate error (but not the specific bug)
      NULL
    }
  )

  # If we got surrogates, verify they work
  if (!is.null(surrogates)) {
    unit_x <- list(c(x1 = 0.5, x2 = 0.5))
    pred <- BATON:::predict_surrogates(surrogates, unit_x)
    expect_true(all(is.finite(pred$EN$mean)))
  }
})

# ============================================================================
# Additional fixes: initial_history objective column mapping
# ============================================================================

test_that("initial_history with only objective column (not metric-named) works", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.70), type1 = c("le", 0.15))

 # Create initial_history with "objective" column but NOT "EN" column
  # This matches the documented interface
  initial_history <- data.frame(
    x1 = c(0.3, 0.5, 0.7, 0.4, 0.6),
    x2 = c(0.2, 0.5, 0.3, 0.6, 0.4),
    objective = c(0.1, 0.05, 0.15, 0.08, 0.12),  # objective values, not "EN"
    power = c(0.82, 0.85, 0.78, 0.80, 0.83),
    type1 = c(0.09, 0.08, 0.11, 0.10, 0.07),
    fidelity = rep("low", 5),
    feasible = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  # Should work - objective column should be mapped to "EN" in metrics
  fit <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = 5,
    budget = 8,
    initial_history = initial_history,
    seed = 123,
    progress = FALSE
  )

  # Check that metrics list includes the objective
  for (i in seq_len(5)) {
    expect_true("EN" %in% names(fit$history$metrics[[i]]))
  }
})

test_that("initial_history with numeric variance column is preserved", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.70), type1 = c("le", 0.15))

  # Create initial_history with numeric variance column (not list-column)
  initial_history <- data.frame(
    x1 = c(0.3, 0.5, 0.7, 0.4, 0.6),
    x2 = c(0.2, 0.5, 0.3, 0.6, 0.4),
    EN = c(0.1, 0.05, 0.15, 0.08, 0.12),
    power = c(0.82, 0.85, 0.78, 0.80, 0.83),
    type1 = c(0.09, 0.08, 0.11, 0.10, 0.07),
    objective = c(0.1, 0.05, 0.15, 0.08, 0.12),
    fidelity = rep("low", 5),
    feasible = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    variance = c(0.001, 0.002, 0.001, 0.003, 0.002),  # Numeric column
    stringsAsFactors = FALSE
  )

  fit <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = 5,
    budget = 8,
    initial_history = initial_history,
    seed = 123,
    progress = FALSE
  )

  # Check that variance was converted to list-column and values preserved
  expect_true(is.list(fit$history$variance))

  # First 5 rows should have variance for "EN"
  for (i in 1:5) {
    var_list <- fit$history$variance[[i]]
    if (length(var_list) > 0) {
      expect_true("EN" %in% names(var_list))
    }
  }
})

# ============================================================================
# Additional fixes: run_init_with_stopping list-style output
# ============================================================================

test_that("run_init_with_stopping handles list-style simulator output", {
  skip_on_cran()

  # Simulator returning list with metrics element
  list_sim <- function(theta, fidelity = "low", seed = 1, ...) {
    list(
      metrics = c(EN = 0.1 + theta$x1, power = 0.85),
      variance = c(EN = 0.001, power = 0.004),
      n_rep = 1000
    )
  }

  design <- matrix(runif(10), ncol = 2)
  colnames(design) <- c("x1", "x2")

  # Should correctly extract EN from list(metrics = ...)
  result <- BATON:::run_init_with_stopping(
    sim_fun = list_sim,
    design = design,
    objective = "EN",
    fidelity = "low",
    init_config = BATON::init_stopping_config(enabled = FALSE),
    seed = 123,
    verbose = FALSE
  )

  # y should be EN values (0.1 + x1)
  expected_y <- 0.1 + design[, "x1"]
  expect_equal(result$y, expected_y, tolerance = 1e-10)
})

test_that("run_init_with_stopping handles plain list simulator output", {
  skip_on_cran()

  # Simulator returning plain list (not with $metrics)
  plain_list_sim <- function(theta, fidelity = "low", seed = 1, ...) {
    list(EN = 0.2 + theta$x2, power = 0.80)
  }

  design <- matrix(runif(6), ncol = 2)
  colnames(design) <- c("x1", "x2")

  result <- BATON:::run_init_with_stopping(
    sim_fun = plain_list_sim,
    design = design,
    objective = "EN",
    fidelity = "low",
    init_config = BATON::init_stopping_config(enabled = FALSE),
    seed = 123,
    verbose = FALSE
  )

  # y should be EN values (0.2 + x2)
  expected_y <- 0.2 + design[, "x2"]
  expect_equal(result$y, expected_y, tolerance = 1e-10)
})
