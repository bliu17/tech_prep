#tags: time-series, regression, stationarity, cointegration, OLS

# Spurious Regression

## Definition

**Spurious regression** occurs when two non-stationary time series appear statistically related — high R², significant t-statistics — purely because they share a common trend, not because of any true causal or economic relationship.

Standard OLS inference assumes stationarity. With integrated (I(1)) series, the regression residuals are also non-stationary, violating OLS assumptions. As a result, t-stats and F-stats do not follow their usual distributions, rendering p-values meaningless.

## Classic Example

Regressing two independent random walks on each other. Despite no true relationship, OLS typically produces a "significant" coefficient and high R².

## How to Detect and Avoid

1. **Test for stationarity** — ADF test, KPSS test
2. **If both series are I(1), test for cointegration**
   - Engle-Granger two-step test
   - Johansen test
   - If cointegrated → regression is valid, residuals are I(0)
   - If not cointegrated → regression is spurious
3. **If not cointegrated, difference the series** before regressing

## Cointegration

Two I(1) series are **cointegrated** if a linear combination of them is I(0) — i.e., they share a common stochastic trend and cannot drift apart permanently. The residual from regressing one on the other is stationary, meaning the relationship is genuine and inference is valid.

**Examples:**

- **Stock price and its futures price** — both are I(1), but the basis (spot − futures) is mean-reverting. Arbitrage forces them to move together.
- **10-year and 2-year Treasury yields** — both wander as I(1) processes, but the spread (10Y − 2Y) is stationary. Monetary policy and market expectations tie them to a common long-run level.

## Related Notes

- [[linear_regression]]
- [[adf_test]]
