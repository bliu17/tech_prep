#tags: #regression #GLS #robust-SE #newey-west #white #sandwich #inference

# GLS, and Its Robust Approximations (White, Newey-West, Fama-MacBeth)

This note collects the Generalized Least Squares framework and the robust-SE methods that approximate it. All of them are answers to the same question: *what do you do when the errors aren't spherical, i.e. $\text{Var}(\varepsilon \mid X) \neq \sigma^2 I$?* See [[linear_regression]] §4 (homoskedasticity) and §5 (autocorrelation) for where these violations arise.

## The Sandwich Variance

For **any** error covariance $\Omega = \text{Var}(\varepsilon \mid X)$, the exact variance of OLS $\hat\beta$ is the sandwich form:

$$\text{Var}(\hat{\beta} \mid X) = (X^\top X)^{-1} \, \underbrace{X^\top \Omega X}_{\text{meat}} \, (X^\top X)^{-1}.$$

The "bread" $(X^\top X)^{-1}$ is fixed; all the action is in estimating the $p \times p$ **meat** $X^\top \Omega X$. White and Newey-West differ only in how they approximate it. None of them ever form or invert the full $\Omega$.

## White (HC0)

Assumes $\Omega$ is **diagonal** (heteroskedasticity only). Approximates each diagonal entry $\sigma_i^2$ with the squared OLS residual $\hat{\varepsilon}_i^2$:

$$X^\top \hat\Omega X = \sum_{i=1}^n \hat{\varepsilon}_i^2 \, x_i x_i^\top.$$

A single $\hat{\varepsilon}_i^2$ is a terrible point estimate of $\sigma_i^2$, but the *sum* converges to $\sum_i \sigma_i^2 x_i x_i^\top$ by a law of large numbers — so HC0 is **consistent**, not unbiased. (Finite-sample variants HC1–HC3 rescale the residuals; see [[linear_regression]] §4.4.)

## Newey-West (HAC)

Same sandwich $(X^\top X)^{-1}(X^\top\hat\Omega X)(X^\top X)^{-1}$ as White, but now $\Omega$ has nonzero off-diagonals (autocorrelation). It approximates the **diagonal** entries with $\hat\varepsilon_t^2$ (exactly White's trick) and the **off-diagonal** entries with the lagged residual product $\hat\varepsilon_t\hat\varepsilon_{t-l} \approx \text{Cov}(\varepsilon_t,\varepsilon_{t-l})$, summed up to lag $L$ with Bartlett weights $w_l = 1 - \tfrac{l}{L+1}$ (which keep the estimate positive semidefinite):

$$X^\top\hat\Omega X = \sum_t \hat\varepsilon_t^2 x_t x_t^\top + \sum_{l=1}^{L} w_l \sum_t \hat\varepsilon_t\hat\varepsilon_{t-l}\big(x_t x_{t-l}^\top + x_{t-l}x_t^\top\big).$$

So White is the $L=0$ case. Like White, it fixes inference only ($\hat\beta$ stays OLS), and never inverts $\Omega$.

## GLS — the General Solution

Given the error covariance $\text{Var}(\varepsilon \mid X) = \Omega$:

$$\hat{\beta}_{GLS} = (X^\top \Omega^{-1} X)^{-1} X^\top \Omega^{-1} y$$

This is BLUE when $\Omega$ is known, recovering **efficiency** (not just valid SEs). WLS, White, and Fama-MacBeth are all special cases or feasible approximations of GLS: **WLS** = GLS with diagonal $\Omega = \text{diag}(\sigma_i^2)$; **White** approximates only GLS's *inference*, plugging $\hat\Omega$ into the sandwich meat (non-inverted) while keeping OLS $\hat\beta$; **Fama-MacBeth** sidesteps the infeasible $N\times N$ cross-sectional $\Omega$ by aggregating in the time dimension.

## The Unifying Picture

| Method | $\Omega$ structure | Uses $\Omega$ how | Fixes |
|---|---|---|---|
| OLS (classical) | $\sigma^2 I$ | — | — |
| White / HC | diagonal | meat, non-inverted | inference |
| Newey-West / HAC | diagonal + banded off-diagonals | meat, non-inverted | inference |
| Fama-MacBeth | block-diagonal (cross-section) | sidestepped via time aggregation | inference |
| **GLS** | full | **inverted** as weights | **estimator (efficiency)** |

The dividing line: robust methods (White, NW, FM) use $\Omega$ **non-inverted** and only repair standard errors, so they tolerate a noisy $\hat\Omega$. GLS uses $\Omega^{-1}$ as weights to repair $\hat\beta$ itself — efficient when $\Omega$ is known, but a noisy estimated $\hat\Omega^{-1}$ (FGLS) propagates violently, which is why robust methods are preferred when the covariance structure isn't credibly known.

## Related Notes

- [[linear_regression]] — §4.4 (White/HC0), §5.4 (Newey-West, reference to this note)
- [[fama_macbeth]]
- [[adf_test]] — HAC assumes errors are stationary / short-memory; non-stationary errors need differencing first
