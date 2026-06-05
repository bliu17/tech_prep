#tags: CLT, central-limit-theorem, asymptotics, sampling-distribution, LLN, normal, statistics

# Central Limit Theorem (CLT)

## Statement

Let $X_1, X_2, \dots, X_n$ be i.i.d. with finite mean $\mu$ and finite variance $\sigma^2 < \infty$. Then the standardized sample mean converges **in distribution** to a standard normal:

$$\frac{\bar{X}_n - \mu}{\sigma / \sqrt{n}} \;\xrightarrow{d}\; N(0, 1) \qquad \text{as } n \to \infty$$

Equivalently, for large $n$:

$$\bar{X}_n \;\approx\; N\!\left(\mu, \frac{\sigma^2}{n}\right), \qquad \sum_{i=1}^n X_i \;\approx\; N\!\left(n\mu, \, n\sigma^2\right)$$

The remarkable part: this holds **regardless of the shape of the underlying distribution** (skewed, discrete, bounded — anything with finite variance). This universality is why the normal distribution is everywhere in inference.

---

## The $\sqrt{n}$ scaling (the key intuition)

The sample mean's *variance* shrinks like $\sigma^2/n$, so its *standard deviation* shrinks like $\sigma/\sqrt{n}$. To get a stable, non-degenerate limit you must rescale by $\sqrt{n}$ — multiply the deviation $(\bar{X}_n - \mu)$ by $\sqrt{n}$:

$$\sqrt{n}\,(\bar{X}_n - \mu) \;\xrightarrow{d}\; N(0, \sigma^2)$$

Too little scaling ($\bar{X}_n - \mu$ alone) collapses to a point mass at 0 (that's the LLN); the right $\sqrt{n}$ scaling reveals the Gaussian fluctuations around the mean. The standard error $\sigma/\sqrt{n}$ falling at rate $1/\sqrt{n}$ — not $1/n$ — is *the* reason quadrupling your sample only halves your error.

---

## LLN vs CLT

These are complementary, often confused:

| | Law of Large Numbers | Central Limit Theorem |
|---|---|---|
| **What it says** | $\bar{X}_n \to \mu$ (converges to a constant) | $\sqrt{n}(\bar{X}_n - \mu) \to N(0,\sigma^2)$ (shape of fluctuations) |
| **Type** | convergence in probability / a.s. | convergence in distribution |
| **Tells you** | *where* the mean lands | *how it's distributed around* that point |

LLN says the average settles down; CLT describes the residual wobble at scale $\sqrt{n}$.

---

## Why it matters (uses across the vault)

The CLT is the engine behind almost all large-sample inference:

- **z-test / t-test** — the sample mean is approximately normal even when the data isn't, justifying these tests for large $n$. See [[hypothesis_testing]] and [[distribution]].
- **OLS inference** — $\hat\beta$ is asymptotically normal, so t- and F-stats are approximately valid even when errors are non-normal (assumption 6 is "not strictly needed"). See [[linear_regression]].
- **MLE asymptotics** — the score is a sum of i.i.d. mean-zero terms; CLT gives the asymptotic normality of $\hat\theta$. See [[maximum_likelihood]].
- **Fama-MacBeth** — the cross-sectional coefficients are averaged over time and their SE comes straight from the CLT. See [[fama_macbeth]].

---

## Common variants

**CLT for a proportion (Bernoulli special case).** If $\hat{p}$ is the sample proportion from $n$ Bernoulli($p$) trials, then since $\text{Var}(X)=p(1-p)$:

$$\hat{p} \;\approx\; N\!\left(p, \frac{p(1-p)}{n}\right)$$

This is the basis of proportion z-tests and confidence intervals. Rule of thumb: the normal approximation is decent when $np \ge 10$ and $n(1-p) \ge 10$.

**Delta method.** For a smooth function $g$ of an asymptotically normal estimator, a first-order Taylor expansion propagates the normality:

$$\sqrt{n}\big(g(\bar{X}_n) - g(\mu)\big) \;\xrightarrow{d}\; N\!\big(0,\; g'(\mu)^2 \sigma^2\big)$$

Used constantly to get standard errors of transformed estimates (ratios, log-returns, Sharpe ratios).

**Multivariate CLT.** For random vectors, $\sqrt{n}(\bar{\mathbf{X}}_n - \boldsymbol\mu) \xrightarrow{d} N(\mathbf{0}, \Sigma)$ with covariance matrix $\Sigma$ — underpins joint inference and the Wald/F frameworks.

**Lindeberg-Lévy** is the classic i.i.d. version above. **Lyapunov / Lindeberg CLTs** relax i.i.d. to merely independent (not identically distributed) under a condition that no single term dominates the variance.

---

## Rate of convergence & when it fails

**Berry-Esseen theorem** bounds the approximation error: the max gap between the true CDF of the standardized mean and the normal CDF is

$$\le \frac{C\,\rho}{\sigma^3 \sqrt{n}}, \qquad \rho = E|X-\mu|^3$$

So convergence is at rate $1/\sqrt{n}$, and **heavier skew/kurtosis (larger third moment $\rho$) means slower convergence** — you need a bigger $n$ before normality kicks in.

**When the CLT does NOT apply:**

- **Infinite / undefined variance.** Fat-tailed distributions like the **Cauchy** (no finite mean or variance) violate the premise — the sample mean of Cauchy variables is *itself* Cauchy, never normal. The relevant generalization is the **Generalized CLT**, whose limits are stable (Lévy) distributions.
- **Strong dependence.** The i.i.d. CLT fails under heavy autocorrelation; you need CLTs for dependent data (mixing conditions), which is exactly why time-series inference uses HAC/Newey-West-style corrections.
- **Finance caveat.** Asset returns have **fat tails and volatility clustering**, so the variance may be large/unstable and dependence is real. Sample means can look "normal enough" while *extreme events are badly underestimated* by a naive CLT argument — a core reason risk models lean on Student-$t$ innovations and robust/heavy-tailed methods rather than assuming Gaussianity.

---

## Related Notes

- [[distribution]]
- [[hypothesis_testing]]
- [[maximum_likelihood]]
- [[linear_regression]]
- [[fama_macbeth]]
- [[variance_covariance_correlation]]
