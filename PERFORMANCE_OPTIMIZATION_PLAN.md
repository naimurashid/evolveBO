# evolveBO Performance Optimization Plan

Based on profiling analysis (1.44s total runtime, lightweight test):
- **76-77%** in `predict_surrogates` / `acq_eci` (GP prediction)
- **50-64%** in `as.data.frame` operations (data coercion)
- **12.5%** in `deparse` / string manipulation

## Priority 1: CRITICAL - `predict_surrogates()` Optimization

**File:** `R/surrogates.R` (lines 424-521)

### Problem
Per-candidate data frame construction with rbind for 2000+ candidates, repeated for each metric:

```r
# Current (SLOW) - lines 448-476
design_list <- lapply(seq_along(unit_x), function(i) {
  point <- unit_x[[i]]
  as.data.frame(as.list(vec[param_names]), stringsAsFactors = FALSE)
})
design <- do.call(rbind, design_list)
```

### Solution
Build design matrix once as a matrix, convert to data.frame only once:

```r
# Optimized - replace lines 448-476
predict_surrogates <- function(surrogates, unit_x, ...) {
  if (length(unit_x) == 0) return(list())

  param_names <- names(unit_x[[1]])
  n_candidates <- length(unit_x)
  n_params <- length(param_names)

  # Build matrix directly (avoid per-candidate data.frame)
  design_mat <- matrix(NA_real_, nrow = n_candidates, ncol = n_params)
  colnames(design_mat) <- param_names

  for (i in seq_len(n_candidates)) {
    design_mat[i, ] <- unlist(unit_x[[i]][param_names])
  }

  # Convert to data.frame ONCE
  design <- as.data.frame(design_mat, stringsAsFactors = FALSE)

  # Predict all metrics using the SAME design matrix
  predictions <- lapply(names(surrogates), function(metric) {
    model <- surrogates[[metric]]
    if (inherits(model, "km")) {
      pred <- predict(model, newdata = design, type = "UK", checkNames = FALSE)
      list(mean = pred$mean, sd = pred$sd)
    } else if (inherits(model, "hetGP")) {
      pred <- predict(model, x = as.matrix(design))
      list(mean = pred$mean, sd = sqrt(pred$sd2))
    } else {
      # Constant predictor fallback
      list(mean = rep(model$mean, n_candidates),
           sd = rep(model$sd, n_candidates))
    }
  })
  names(predictions) <- names(surrogates)
  predictions
}
```

**Expected improvement:** 40-60% reduction in `predict_surrogates` time

---

## Priority 2: HIGH - `record_evaluation()` Optimization

**File:** `R/bo_calibrate.R` (lines 745-869)

### Problem 2a: Metrics unpacking (lines 786-793)
```r
# Current (SLOW)
metrics_list <- lapply(as.list(metrics), function(x) {
  if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)
})
as.data.frame(metrics_list, stringsAsFactors = FALSE)
```

### Solution 2a: Direct vector to data.frame
```r
# Optimized
metrics_vec <- vapply(metrics, function(x) {
  if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)
}, numeric(1))
metrics_df <- as.data.frame(as.list(metrics_vec), stringsAsFactors = FALSE)
```

### Problem 2b: Theta unpacking (lines 817-823)
Same pattern as metrics - replace with vapply.

### Problem 2c: Growing history with bind_rows (line 853)
```r
# Current (SLOW) - O(n^2) growth
history <- dplyr::bind_rows(history, new_row)
```

### Solution 2c: Pre-allocate history, fill in place
```r
# At initialization (around line 280):
max_evals <- budget + n_init + 10  # buffer
history <- data.frame(
  eval_id = rep(NA_integer_, max_evals),
  iteration = rep(NA_integer_, max_evals),
  # ... pre-allocate all columns
)
history_idx <- 0L

# In record_evaluation:
history_idx <<- history_idx + 1L
history[history_idx, ] <- new_row_values
# At end, trim: history <- history[seq_len(history_idx), ]
```

**Expected improvement:** 20-30% reduction in data frame overhead

---

## Priority 3: MEDIUM - `lhs_candidate_pool()` Optimization

**File:** `R/bo_calibrate.R` (lines 974-982)

### Problem
Creates list of named vectors for 2000+ candidates:
```r
# Current (SLOW)
lapply(seq_len(n), function(i) {
  values <- lhs[i, , drop = TRUE]
  names(values) <- names(bounds)
  values
})
```

### Solution
Return matrix with column names, convert to list only when needed:
```r
# Optimized
lhs_candidate_pool <- function(n, bounds) {
  lhs_mat <- lhs::randomLHS(n, length(bounds))
  colnames(lhs_mat) <- names(bounds)
  # Return matrix for internal use
  lhs_mat
}

# Helper to convert single row to named vector when needed
row_to_named_vec <- function(mat, i) {
  vec <- mat[i, , drop = TRUE]
  names(vec) <- colnames(mat)
  vec
}
```

Then update callers to work with matrix representation internally.

**Expected improvement:** 10-15% reduction in candidate generation overhead

---

## Priority 4: LOW - String Manipulation Cleanup

**File:** `R/bo_calibrate.R` (lines 843-844)

### Problem
Column renaming on every evaluation:
```r
names(metrics_df)[names(metrics_df) %in% conflict_cols] <-
  paste0("metric_", names(metrics_df)[names(metrics_df) %in% conflict_cols])
```

### Solution
Pre-compute renamed column names once at initialization:
```r
# At initialization:
metric_col_map <- setNames(
  paste0("metric_", conflict_cols),
  conflict_cols
)

# In record_evaluation:
if (any(names(metrics_df) %in% names(metric_col_map))) {
  names(metrics_df) <- ifelse(
    names(metrics_df) %in% names(metric_col_map),
    metric_col_map[names(metrics_df)],
    names(metrics_df)
  )
}
```

**Expected improvement:** 5-10% reduction in string operations

---

## Implementation Order

1. **Phase 1 (CRITICAL):** `predict_surrogates()` matrix optimization
   - Highest impact, isolated change
   - Test: Re-run profiler, expect 40-60% improvement in GP prediction

2. **Phase 2 (HIGH):** `record_evaluation()` pre-allocation
   - Moderate complexity, affects history management
   - Test: Verify history integrity after runs

3. **Phase 3 (MEDIUM):** Candidate pool matrix representation
   - Requires updating multiple callers
   - Test: Verify acquisition function correctness

4. **Phase 4 (LOW):** String operation cleanup
   - Minor impact, easy to implement
   - Test: Verify column names are correct

---

## Testing Strategy

1. **Correctness tests:**
   - Run existing test suite after each phase
   - Compare BO results (objective values, selected points) before/after
   - Verify feasibility constraints are correctly evaluated

2. **Performance tests:**
   - Re-run `scripts/profile_bottlenecks_light.R` after each phase
   - Target: <0.5s for lightweight profile (down from 1.44s)
   - Full production run: target 30-50% wall-time reduction

3. **Regression tests:**
   - Run `test_integration.R` with fixed seed
   - Compare final objective values (should be identical with same seed)

---

## Expected Overall Improvement

| Phase | Target Reduction | Cumulative |
|-------|------------------|------------|
| Phase 1 (predict_surrogates) | 40-60% of GP time | ~50% overall |
| Phase 2 (record_evaluation) | 20-30% of data ops | ~60% overall |
| Phase 3 (candidate pool) | 10-15% of generation | ~65% overall |
| Phase 4 (string ops) | 5-10% of misc | ~70% overall |

**Target:** 60-70% reduction in BO overhead (non-simulator time)
