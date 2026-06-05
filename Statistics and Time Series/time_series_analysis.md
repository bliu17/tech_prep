#tags: time-series, AR, MA, ARMA, ARIMA, ACF, PACF, ARCH, GARCH, stationarity, volatility

# Time Series Analysis

## Contents

- [[#Stationarity]]
- [[#White Noise]]
- [[#AR (Autoregressive)]]
- [[#MA (Moving Average)]]
- [[#ARMA]]
- [[#ARIMA]]
- [[#SARIMA (Seasonal)]]
- [[#ACF and PACF — Identification Toolkit]]
- [[#Box-Jenkins Model-Building Workflow]]
- [[#ARCH — Modeling Volatility]]
- [[#GARCH]]
- [[#Identifying ARCH/GARCH Orders]]

---

## Stationarity

A series is **(weakly / covariance) stationary** if its first two moments are time-invariant:

$$E[y_t] = \mu, \qquad \text{Var}(y_t) = \gamma_0, \qquad \text{Cov}(y_t, y_{t-k}) = \gamma_k \;\text{(depends only on lag } k)$$

Stationarity is the prerequisite for ARMA modeling — coefficients are only meaningful if the process's statistical properties don't drift. Test with the [[adf_test]] (ADF, null = unit root) and KPSS (null = stationary); run both. Non-stationary series are made stationary by **differencing** (for stochastic trends / unit roots) or **detrending** (for deterministic trends). Over-differencing is a real cost — it injects an artificial MA unit root and inflates variance.

---

## White Noise

The building block: $\varepsilon_t \sim \text{WN}(0, \sigma^2)$ with $E[\varepsilon_t]=0$, constant variance, and zero autocorrelation at all lags. ACF and PACF are both **flat at zero** (within confidence bands). A well-specified model leaves white-noise residuals — that's the diagnostic target.

---

## AR (Autoregressive)

Today's value is a regression on its own past:

$$\text{AR}(p): \quad y_t = c + \phi_1 y_{t-1} + \phi_2 y_{t-2} + \cdots + \phi_p y_{t-p} + \varepsilon_t$$

**Stationarity condition:** the roots of the characteristic polynomial $1 - \phi_1 z - \cdots - \phi_p z^p = 0$ must lie **outside** the unit circle. For AR(1) this reduces to $|\phi_1| < 1$.

**AR(1) moments:** $E[y_t] = \frac{c}{1-\phi}$, $\text{Var}(y_t) = \frac{\sigma^2}{1-\phi^2}$, and $\rho_k = \phi^k$ — the ACF **decays geometrically**.

**Signature:** ACF tails off (decays), PACF cuts off after lag $p$.

---

## MA (Moving Average)

Today's value is a weighted sum of current and past *shocks*:

$$\text{MA}(q): \quad y_t = \mu + \varepsilon_t + \theta_1 \varepsilon_{t-1} + \cdots + \theta_q \varepsilon_{t-q}$$

**Always stationary** (a finite sum of white-noise terms has constant moments), so no stationarity restriction on $\theta$. For meaningful estimation you instead want **invertibility** (roots of the MA polynomial outside the unit circle) so the process has a unique AR($\infty$) representation.

**MA(1) facts:** autocorrelation is **zero beyond lag 1** — $\rho_1 = \frac{\theta_1}{1+\theta_1^2}$, $\rho_k = 0$ for $k > 1$. Note $|\rho_1| \le 0.5$ always for MA(1).

**Signature:** ACF cuts off after lag $q$, PACF tails off (decays).

---

## ARMA

Combines both — captures series with autoregressive memory *and* shock dependence:

$$\text{ARMA}(p,q): \quad y_t = c + \sum_{i=1}^p \phi_i y_{t-i} + \varepsilon_t + \sum_{j=1}^q \theta_j \varepsilon_{t-j}$$

ARMA is parsimonious: a high-order pure-AR or pure-MA model can often be replaced by a low-order ARMA. **Signature:** *both* ACF and PACF tail off — neither cuts cleanly, which is exactly why mixed models are harder to identify by eye (see EACF / information criteria below).

---

## ARIMA

ARIMA($p,d,q$) = ARMA($p,q$) applied to a series differenced $d$ times. The "I" is **Integrated**: $d$ is the number of differences needed to reach stationarity.

$$\nabla^d y_t = (1-L)^d y_t \;\;\text{is ARMA}(p,q), \qquad L = \text{lag operator } (L y_t = y_{t-1})$$

- $d=0$ → already stationary (plain ARMA)
- $d=1$ → one difference removes a single unit root / linear trend (most common for prices, levels)
- $d=2$ → rarely needed; for series with a changing trend

Determine $d$ first (via ADF/KPSS), *then* identify $p,q$ on the differenced series.

---

## SARIMA (Seasonal)

SARIMA$(p,d,q)(P,D,Q)_s$ adds a seasonal block with period $s$ (e.g. $s=12$ monthly, $s=4$ quarterly). The seasonal $(P,D,Q)$ terms act on lags that are multiples of $s$. Diagnostic tell: **spikes in the ACF at seasonal lags** ($s, 2s, 3s, \dots$). Seasonal differencing $(1-L^s)$ handles seasonal unit roots.

---

## ACF and PACF — Identification Toolkit

This is the core interview skill: read the two correlograms to pick $p$ and $q$.

**ACF (Autocorrelation Function):** correlation between $y_t$ and $y_{t-k}$ — includes *all* intermediate effects (the "total" correlation at lag $k$).

**PACF (Partial Autocorrelation Function):** correlation between $y_t$ and $y_{t-k}$ **after removing** the linear effect of the intermediate lags $y_{t-1}, \dots, y_{t-k+1}$ — the "direct" effect. Conceptually it's the last coefficient $\phi_{kk}$ from regressing $y_t$ on its first $k$ lags (same partialling-out idea as [[Frisch-Waugh-Lovell]]).

### The identification table

| Model | ACF | PACF |
|---|---|---|
| **AR($p$)** | tails off (decays) | **cuts off after lag $p$** |
| **MA($q$)** | **cuts off after lag $q$** | tails off (decays) |
| **ARMA($p,q$)** | tails off | tails off |
| White noise | no spikes | no spikes |

**Mnemonic:** *AR → cut in PACF; MA → cut in ACF.* The cut tells you the order; the cutting plot tells you the type. (AR order is read off PACF, MA order off ACF — they cross.)

**Why the cutoffs happen:**
- An AR($p$) has only $p$ direct lag dependencies, so its PACF is exactly zero beyond lag $p$.
- An MA($q$) is a finite window of $q+1$ shocks, so observations more than $q$ apart share no shocks → ACF is zero beyond $q$.

**Reading the bands:** a spike is "significant" if it pierces the $\pm 1.96/\sqrt{n}$ confidence band. Decay can be geometric or damped-sinusoidal (complex AR roots).

---

## Box-Jenkins Model-Building Workflow

1. **Stationarize** — plot, ADF/KPSS, difference until stationary → fixes $d$.
2. **Identify** — inspect ACF/PACF to propose $p, q$ (and seasonal $P,Q$).
3. **Estimate** — fit by MLE / conditional least squares (see [[maximum_likelihood]]).
4. **Select** — compare candidates by **AIC** ($-2\ln L + 2k$) or **BIC** ($-2\ln L + k\ln n$); BIC penalizes complexity harder → more parsimonious models. For mixed ARMA where ACF/PACF are ambiguous, the **EACF (Extended ACF)** or a grid search over $(p,q)$ minimizing AIC/BIC is the practical tool.
5. **Diagnose** — residuals should be white noise: check residual ACF and the **Ljung-Box Q-test** ($H_0$: residuals jointly uncorrelated). If it fails, revise.
6. **Forecast.**

---

## ARCH — Modeling Volatility

ARMA models the *mean*; ARCH/GARCH model the **conditional variance**. Financial returns are typically near-white-noise in the mean but exhibit **volatility clustering** — big moves follow big moves (Mandelbrot). The series is uncorrelated but *not independent*: $r_t$ has no autocorrelation, yet $r_t^2$ does.

**ARCH($q$)** — the conditional variance depends on past squared shocks:

$$r_t = \mu + \varepsilon_t, \quad \varepsilon_t = \sigma_t z_t, \;\; z_t \sim \text{iid}(0,1)$$
$$\sigma_t^2 = \omega + \alpha_1 \varepsilon_{t-1}^2 + \cdots + \alpha_q \varepsilon_{t-q}^2$$

Constraints for a valid (positive) variance: $\omega > 0$, $\alpha_i \ge 0$; for stationarity $\sum \alpha_i < 1$. ARCH reproduces clustering and **fat tails** (excess kurtosis) even with Gaussian $z_t$. Drawback: often needs a large $q$ to fit the slow decay of volatility memory.

---

## GARCH

**GARCH($p,q$)** (Bollerslev) generalizes ARCH by letting variance also depend on its own past values — far more parsimonious:

$$\sigma_t^2 = \omega + \sum_{i=1}^q \alpha_i \varepsilon_{t-i}^2 + \sum_{j=1}^p \beta_j \sigma_{t-j}^2$$

- $\alpha$ (ARCH term) = reaction to recent shocks (short-run jumps in vol).
- $\beta$ (GARCH term) = persistence / memory of volatility.
- **Constraints:** $\omega>0$, $\alpha_i, \beta_j \ge 0$; **stationarity / mean-reversion:** $\sum \alpha_i + \sum \beta_j < 1$.
- **Unconditional variance:** $\bar{\sigma}^2 = \dfrac{\omega}{1 - \sum\alpha_i - \sum\beta_j}$. The closer $\sum\alpha+\sum\beta$ is to 1, the more persistent (slowly mean-reverting) the volatility. If it equals 1, you have **IGARCH** (integrated, infinite memory).

**GARCH(1,1) is the workhorse** — fits most financial return series well:
$$\sigma_t^2 = \omega + \alpha_1 \varepsilon_{t-1}^2 + \beta_1 \sigma_{t-1}^2$$

**Common extensions** (for the leverage effect — negative returns raise vol more than positive ones, which symmetric GARCH misses): **EGARCH** (models $\ln\sigma_t^2$, no positivity constraints needed) and **GJR-GARCH / TGARCH** (adds an asymmetric term that switches on for negative shocks).

---

## Identifying ARCH/GARCH Orders

Volatility orders are **not** read off the ordinary ACF/PACF of the series (those govern the mean). The procedure:

1. **Fit the mean model first** (constant, or ARMA) and take residuals $\hat{\varepsilon}_t$.
2. **Test for ARCH effects** on the residuals:
   - Plot the **ACF/PACF of squared residuals $\hat{\varepsilon}_t^2$** — significant autocorrelation signals ARCH/GARCH is needed.
   - **Engle's ARCH-LM test:** regress $\hat{\varepsilon}_t^2$ on $q$ of its own lags; $n R^2 \sim \chi^2_q$ under $H_0$ (no ARCH). Reject → volatility clustering present.
3. **Choose orders:** the ACF/PACF of $\hat{\varepsilon}_t^2$ suggest the ARCH order $q$ analogously to the mean case, but in practice **GARCH(1,1) is the default starting point** and orders are selected by **AIC/BIC**. Higher orders rarely help — most return series are adequately captured by GARCH(1,1).
4. **Estimate** by maximum likelihood (see [[maximum_likelihood]]) — typically assuming Gaussian or Student-$t$ innovations ($t$ captures the remaining fat tails better).
5. **Re-check:** the **standardized** residuals $\hat{\varepsilon}_t/\hat{\sigma}_t$ and their squares should now be white noise (re-run ARCH-LM / Ljung-Box).

**One-line summary:** use ACF/PACF of the *level* to pick ARMA orders for the mean; use ACF/PACF of *squared residuals* (plus ARCH-LM) to detect and size the GARCH variance model.

---

## Related Notes

- [[adf_test]]
- [[spurious_regression]]
- [[maximum_likelihood]]
- [[distribution]]
- [[GLS]]
