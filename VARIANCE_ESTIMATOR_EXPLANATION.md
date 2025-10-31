# Variance Estimator Explanation

## The Question
Why does the default variance estimator return `NA` for continuous metrics like `EN` (expected sample size) and `ET` (expected trial duration)?

## Background: What is Variance Here?

When you run a clinical trial simulator with, say, 1000 Monte Carlo replications, you get:
- **Mean estimates**: `power = 0.85`, `EN = 245.7`
- **Variance estimates**: How uncertain are these means due to Monte Carlo sampling?

This variance is **not** the variance of the underlying distribution, but rather the **standard error of the Monte Carlo estimate**.

## The Two Cases

### Case 1: Proportion Metrics (values in [0, 1])

**Examples**: `power`, `type1` (Type I error rate)

These are proportions estimated from binary outcomes across replications:
- Power = proportion of simulations where null hypothesis was rejected
- Type I error = proportion of simulations with false positive

**Variance Formula**: For a proportion `p` estimated from `n` replications:
```
Var(p̂) = p(1-p) / n
```

This is the **binomial variance formula**, which is mathematically exact for proportions.

**Example**:
```r
power = 0.85
n_rep = 1000
variance = 0.85 * (1 - 0.85) / 1000 = 0.0001275
standard_error = sqrt(0.0001275) ≈ 0.0113
```

So we estimate power = 0.85 ± 0.011 (95% CI).

---

### Case 2: Continuous Metrics (values outside [0, 1])

**Examples**: `EN` (expected sample size), `ET` (expected trial duration)

These are **continuous random variables**, not proportions:
- `EN = 245.7` might come from averaging actual sample sizes across simulations
- `ET = 18.3` months might come from averaging trial durations

**The Problem**: We cannot estimate variance without knowing the underlying distribution!

For proportions, we know it's binomial, so `p(1-p)/n` works. But for continuous metrics:
- Sample size might follow a complex distribution (depends on adaptive design rules)
- Trial duration might depend on enrollment rates, event times, etc.

**Why we can't use a formula**:
```r
# This doesn't work for continuous metrics:
EN = 245.7
variance = EN * (1 - EN) / n_rep  # ❌ Makes no sense!
# Gives: 245.7 * (-244.7) / 1000 = -60.1 (negative variance?!)
```

---

## The Solution: Return NA and Use a Nugget

When variance is `NA`, the code in `surrogates.R` does this:

```r
if (all(is.na(noise_vec))) {
  nugget <- 1e-6        # Use small nugget
  noise_vec <- NULL     # Don't use heteroskedastic noise
}
```

### What is a "nugget"?

A **nugget** is a small constant added to the diagonal of the GP covariance matrix:

```
Cov(y_i, y_j) = k(x_i, x_j) + δ_ij * nugget
```

Where `δ_ij = 1` if `i=j`, else 0.

**Effect**:
- **Homoskedastic**: Assumes all observations have the same noise level (1e-6)
- **Regularization**: Prevents numerical issues when points are very close together
- **Not optimal, but safe**: Better than crashing or using a wrong variance estimate

---

## What Happens in Each Case?

### Scenario A: Simulator provides variance (BEST)

```r
sim_fun <- function(theta, ...) {
  # Run 1000 simulations
  results <- replicate(1000, run_one_trial(theta))

  res <- c(
    power = mean(results$reject_null),
    EN = mean(results$sample_size)
  )

  # Provide actual Monte Carlo variances
  attr(res, "variance") <- c(
    power = var(results$reject_null) / 1000,
    EN = var(results$sample_size) / 1000
  )
  attr(res, "n_rep") <- 1000

  return(res)
}
```

**Result**:
- Uses **heteroskedastic GP** with actual noise variances
- Optimal uncertainty quantification ✅

---

### Scenario B: Simulator doesn't provide variance (FALLBACK)

```r
sim_fun <- function(theta, ...) {
  # Run 1000 simulations
  results <- replicate(1000, run_one_trial(theta))

  res <- c(
    power = mean(results$reject_null),  # 0.85
    EN = mean(results$sample_size)      # 245.7
  )

  # No variance attribute provided
  return(res)
}
```

**What happens**:

1. `default_variance_estimator()` is called:
   ```r
   power: value = 0.85 (in [0,1]) → variance = 0.85*(1-0.85)/1000 = 0.0001275
   EN: value = 245.7 (not in [0,1]) → variance = NA
   ```

2. GP fitting for `power`:
   ```r
   noise_vec = [0.0001275]  # Has a value
   → Uses heteroskedastic GP with noise.var parameter ✅
   ```

3. GP fitting for `EN`:
   ```r
   noise_vec = [NA]
   all(is.na(noise_vec)) = TRUE
   → nugget = 1e-6
   → noise_vec = NULL
   → Uses homoskedastic GP with small nugget ⚠️
   ```

---

## Why This is Reasonable

1. **For proportions**: We get the mathematically correct variance (binomial formula)

2. **For continuous metrics**:
   - We can't estimate variance without the raw data or knowing the distribution
   - Using a small nugget (1e-6) is a safe default that:
     - Prevents numerical instability
     - Allows GP fitting to succeed
     - Assumes constant (homoskedastic) noise
   - This is **not optimal** but better than:
     - ❌ Crashing
     - ❌ Using a wrong variance formula
     - ❌ Assuming zero noise (overfitting)

---

## Recommendation

**Best practice**: Have your simulator return variance estimates:

```r
sim_fun <- function(theta, fidelity = "high", seed = NULL, ...) {
  n_rep <- switch(fidelity, low = 200, med = 1000, high = 10000)

  set.seed(seed)
  results <- replicate(n_rep, {
    # Simulate one trial
    trial <- run_adaptive_trial(theta)
    c(
      power = trial$reject_null,
      type1 = trial$false_positive,
      EN = trial$sample_size,
      ET = trial$duration_months
    )
  })

  # Compute means
  metrics <- rowMeans(results)

  # Compute Monte Carlo variances (variance of the mean)
  variance <- apply(results, 1, var) / n_rep

  attr(metrics, "variance") <- variance
  attr(metrics, "n_rep") <- n_rep

  return(metrics)
}
```

This gives you:
- Correct variance for ALL metrics (proportions and continuous)
- Heteroskedastic GP that accounts for different noise levels
- Better optimization performance

---

## Summary Table

| Metric | Type | Value | Default Variance Estimate | GP Type |
|--------|------|-------|--------------------------|---------|
| `power` | Proportion | 0.85 | `0.85*(1-0.85)/n` ✅ | Heteroskedastic |
| `type1` | Proportion | 0.05 | `0.05*(1-0.05)/n` ✅ | Heteroskedastic |
| `EN` | Continuous | 245.7 | `NA` → nugget ⚠️ | Homoskedastic |
| `ET` | Continuous | 18.3 | `NA` → nugget ⚠️ | Homoskedastic |

**Legend**:
- ✅ Mathematically correct variance estimate
- ⚠️ Safe fallback (homoskedastic with nugget)

---

## Conclusion

The `NA` return for continuous metrics is **intentional and correct** because:

1. We **cannot** mathematically estimate Monte Carlo variance for arbitrary continuous distributions without raw data
2. Returning `NA` signals "I don't know the variance"
3. The GP fitting code gracefully handles this by using a small nugget (homoskedastic noise assumption)
4. This is safer than guessing or crashing

**The fix applied**: Added comprehensive documentation explaining this behavior so users understand they should provide variance attributes for best results.
