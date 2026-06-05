#tags: time-series, stationarity, unit-root, ADF, hypothesis-testing

# Augmented Dickey-Fuller (ADF) Test

## Purpose

Tests whether a time series has a **unit root** — i.e., whether it is non-stationary (I(1)). Used as a prerequisite before regression to avoid spurious results.

## Hypotheses

$$H_0: \text{unit root exists (series is non-stationary, I(1))}$$
$$H_1: \text{no unit root (series is stationary, I(0))}$$

Note the null is non-stationarity — you need strong evidence to reject it.

## The Regression

The ADF test runs the following regression:

$$\Delta y_t = \alpha + \beta t + \gamma y_{t-1} + \sum_{j=1}^{p} \delta_j \Delta y_{t-j} + \varepsilon_t$$

where:
- $\alpha$ = intercept (optional)
- $\beta t$ = deterministic trend (optional)
- $\gamma y_{t-1}$ = the key term — tests for unit root
- $\sum \delta_j \Delta y_{t-j}$ = lagged differences added to remove autocorrelation in residuals (the "augmented" part)

The test statistic is the t-statistic on $\gamma$:

$$H_0: \gamma = 0 \quad \text{(unit root)} \qquad H_1: \gamma < 0 \quad \text{(stationary)}$$

## Critical Values

The test statistic does **not** follow a standard t-distribution under $H_0$ — it follows the Dickey-Fuller distribution, which is left-skewed. Critical values are more negative than standard t-critical values (e.g., ~−2.86 at 5% for the no-trend case).

## Choosing Lag Length $p$

- Too few lags: residual autocorrelation, size distortion
- Too many lags: loss of power
- Common approach: use AIC or BIC to select $p$

## Limitations

- **Low power** in small samples — tends to fail to reject $H_0$ even when the series is stationary
- Sensitive to the choice of deterministic terms (intercept, trend)
- Alternative: **KPSS test**, where $H_0$ is stationarity — useful to run both as a robustness check

## Related Notes

- [[spurious_regression]]
- [[linear_regression]]
