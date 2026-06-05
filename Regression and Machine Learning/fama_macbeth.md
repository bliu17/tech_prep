#tags: #regression #fama-macbeth #panel-data #inference #factor-models

# Fama-MacBeth Regression

## Why Pooled OLS Fails on Panel Data

Suppose we have panel data with $N$ securities, $K$ features, and $T$ time periods. The model is:

$$r_{i,t} = X_{i,t} \beta + \varepsilon_{i,t}$$

If we stack all $N \times T$ observations and run a single pooled OLS, the residuals violate the independence assumption in two ways:

**1. Cross-sectional correlation:** At any time $t$, residuals $\varepsilon_{i,t}$ and $\varepsilon_{j,t}$ are correlated across stocks (due to common factor exposure, macro shocks, etc.). The $N$ observations at time $t$ are not independent draws.

**2. Time-series autocorrelation:** For the same stock $i$, residuals $\varepsilon_{i,t}$ and $\varepsilon_{i,t+s}$ are correlated over time.

Pooled OLS ignores both. As a result:

- The $N \times N$ residual covariance matrix $\Omega$ is misspecified as $\sigma^2 I$
- Standard errors are **underestimated**
- $t$-statistics are **inflated**
- You **over-reject** the null hypothesis

Correcting this properly would require estimating $\Omega$, which has $N(N+1)/2$ free parameters — infeasible for large $N$.

---

## The Fama-MacBeth Solution

Instead of modeling $\Omega$, Fama-MacBeth (1973) sidesteps the problem entirely by **running $T$ separate cross-sectional regressions**, one per time period.

### Step 1: Cross-sectional regressions

At each $t = 1, \ldots, T$, regress returns on features across all $N$ stocks:

$$\hat{\beta}_t = (X_t^\top X_t)^{-1} X_t^\top r_t$$

This gives a $K$-vector of slopes $\hat{\beta}_t$ for each period.

### Step 2: Time-series average

The Fama-MacBeth estimator is the mean of these $T$ estimates:

$$\hat{\beta}_{FM} = \frac{1}{T} \sum_{t=1}^{T} \hat{\beta}_t$$

**Why is this valid?** If $\beta$ is constant over time, each $\hat{\beta}_t$ is an unbiased estimate of the same $\beta$. Averaging $T$ unbiased estimates reduces variance — by the law of large numbers, $\hat{\beta}_{FM} \to \beta$ as $T \to \infty$.

---

## Computing the $t$-Statistic

The key insight is that inference is done **in the time dimension**, not across observations. Treat $\hat{\beta}_1, \hat{\beta}_2, \ldots, \hat{\beta}_T$ as $T$ i.i.d. draws of the same quantity.

### Deriving SE from CLT

Assume $\hat{\beta}_1, \ldots, \hat{\beta}_T$ are i.i.d. with mean $\beta$ and variance $\sigma^2_\beta = \text{Var}(\hat{\beta}_t)$. By [[central_limit_theorem|CLT]]:

$$\sqrt{T}(\hat{\beta}_{FM} - \beta) \xrightarrow{d} \mathcal{N}(0,\ \sigma^2_\beta)$$

Equivalently, $\hat{\beta}_{FM} - \beta \sim \mathcal{N}(0,\ \sigma^2_\beta / T)$, so:

$$\text{Var}(\hat{\beta}_{FM}) = \frac{\sigma^2_\beta}{T} \implies \text{SE}(\hat{\beta}_{FM}) = \frac{\sigma_\beta}{\sqrt{T}}$$

The step is: if $\sqrt{T} \cdot X \sim \mathcal{N}(0, \sigma^2)$, then $X \sim \mathcal{N}(0, \sigma^2/T)$, giving $\text{SD}(X) = \sigma/\sqrt{T}$.

We estimate $\sigma_\beta$ from the sample standard deviation of $\hat{\beta}_t$:

$$\hat{\sigma}_\beta = \sqrt{\frac{1}{T-1} \sum_{t=1}^{T} \left(\hat{\beta}_t - \hat{\beta}_{FM}\right)^2}$$

Plugging in gives the FM standard error:

$$\text{SE}(\hat{\beta}_{FM}) = \frac{\hat{\sigma}_\beta}{\sqrt{T}} = \sqrt{\frac{1}{T(T-1)} \sum_{t=1}^{T} \left(\hat{\beta}_t - \hat{\beta}_{FM}\right)^2}$$

**Key implication:** FM inference is asymptotic in $T$, not $N$. Adding more stocks does not improve inference — adding more time periods does.

### Variance of $\hat{\beta}_{FM}$

$$\widehat{\text{Var}}(\hat{\beta}_{FM}) = \frac{1}{T^2} \sum_{t=1}^{T} \left(\hat{\beta}_t - \hat{\beta}_{FM}\right)^2$$

This is just the sample variance of the time series of slopes, scaled by $1/T$.

### $t$-statistic

This is the $t$-statistic under the null hypothesis $H_0: \beta = 0$. For a general null $H_0: \beta = c$, replace the numerator with $\hat{\beta}_{FM} - c$.

$$t = \frac{\hat{\beta}_{FM}}{\text{SE}(\hat{\beta}_{FM})}$$

This $t$-statistic is valid because the variance is estimated from the **time-series variation** of $\hat{\beta}_t$, which naturally absorbs cross-sectional correlation. Whatever common shock hits all stocks at time $t$ shows up as a large deviation of $\hat{\beta}_t$ from its mean — increasing the standard error appropriately.

---

## What FM Fixes vs. What It Doesn't

| Issue | FM Fixes? | Notes |
|-------|-----------|-------|
| Cross-sectional residual correlation | ✅ Yes | Absorbed into time-series variance of $\hat{\beta}_t$ |
| Time-series autocorrelation in $\hat{\beta}_t$ | ❌ No | If slopes are persistent, SE is still too small |
| Time-varying $\beta$ | ❌ No | FM assumes $\beta$ is constant; average is misleading if not |

### Fix for autocorrelation: Newey-West

If $\hat{\beta}_t$ is autocorrelated (common for slow-moving factors), apply Newey-West to the time series of $\hat{\beta}_t$:

$$\widehat{\text{Var}}_{NW}(\hat{\beta}_{FM}) = \frac{1}{T^2} \left[ \sum_{t} u_t^2 + 2\sum_{l=1}^{L} \left(1 - \frac{l}{L+1}\right) \sum_{t=l+1}^{T} u_t u_{t-l} \right]$$

where $u_t \equiv \hat{\beta}_t - \hat{\beta}_{FM}$ is the demeaned slope at time $t$ (deviation of each period's estimate from the overall mean), and $L$ is the lag truncation. The first term $\sum_t u_t^2$ is the standard variance; the second term adds weighted autocovariances $u_t u_{t-l}$ to correct for autocorrelation in the slope series.

---

## Alternative Methods

### Clustered Standard Errors

Instead of running separate cross-sectional regressions (FM), run one pooled OLS on all $N \times T$ observations and correct the standard errors via a sandwich estimator:

$$\widehat{\text{Var}}_{clustered}(\hat{\beta}) = (X^\top X)^{-1} B (X^\top X)^{-1}$$

The "meat" $B$ allows arbitrary residual correlation within clusters but assumes independence across clusters. The key choice is the clustering dimension:

- **Cluster by time** — allows cross-sectional correlation within each day. Asymptotically equivalent to FM.
- **Cluster by stock** — allows serial correlation within each stock over time.
- **Double cluster (by both)** — handles both simultaneously via inclusion-exclusion: $B = B_{stock} + B_{time} - B_{individual}$.

The score contribution of each observation is $x_{i,t} \hat{u}_{i,t}$ (not just $\hat{u}_{i,t}$), because only residuals paired with large regressors actually influence $\hat{\beta}$.

**Practical caveat:** clustered SEs require enough clusters ($\geq 30$–$50$) in each dimension to be reliable. With few time periods or few stocks, they are poorly estimated.

**FM vs. clustered by time:** both solve the same cross-sectional dependence problem. FM is more transparent (you see period-by-period slopes); clustering is more flexible and integrates naturally into a single regression framework.

## Summary

Fama-MacBeth works by converting a hard panel inference problem (intractable $N \times N$ covariance) into a simple time-series inference problem (variance of $T$ slope estimates). It's the standard approach in empirical asset pricing for testing whether factors predict returns.

## See Also

- [[linear_regression]]
