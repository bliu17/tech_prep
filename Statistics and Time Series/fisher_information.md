#tags: fisher-information, mle, cramer-rao, estimation, statistics

# Fisher Information

How much information a sample carries about a parameter $\theta$ — the curvature of the log-likelihood. Sets the floor on estimator variance (Cramér-Rao) and the asymptotic variance of the MLE.

---

## Definition

For a density $f(x;\theta)$, the **score** is the derivative of the log-likelihood, and Fisher information is the variance of the score:

$$I(\theta) = \mathbb{E}\!\left[\left(\frac{\partial}{\partial\theta}\log f(X;\theta)\right)^2\right] = -\,\mathbb{E}\!\left[\frac{\partial^2}{\partial\theta^2}\log f(X;\theta)\right]$$

The two forms are equal under regularity conditions. The **negative-second-derivative** form is usually the easier one to compute.

For $n$ i.i.d. observations, information **adds**: $I_n(\theta) = n\,I_1(\theta)$.

## Intuition — curvature

Fisher information measures how **sharply peaked** the log-likelihood is around the true $\theta$.

- **High curvature** (sharp peak) → the data strongly distinguishes nearby parameter values → easy to estimate $\theta$ precisely → high information.
- **Flat likelihood** → many $\theta$ values fit the data almost equally well → hard to pin down → low information.

It is a property of **the assumed model + parameter**, not a model-free property of the data. Change the distribution assumption and the information changes.

## Cramér-Rao lower bound

Any **unbiased** estimator $\hat\theta$ satisfies:

$$\text{Var}(\hat\theta) \ge \frac{1}{I_n(\theta)}$$

An estimator that attains this bound is **efficient**. The MLE attains it asymptotically:

$$\hat\theta_{\text{MLE}} \;\approx\; \mathcal{N}\!\left(\theta, \frac{1}{n\,I_1(\theta)}\right)$$

## Worked example — mean return of i.i.d. normal returns

Returns $X_1,\dots,X_n \overset{\text{iid}}{\sim} \mathcal{N}(\mu,\sigma^2)$, with $\sigma$ **known**. Estimate the mean return $\mu$.

**Log-likelihood** (single observation, dropping $\mu$-free constants):

$$\log f(x_i;\mu) = -\frac{(x_i-\mu)^2}{2\sigma^2} + \text{const}$$

**Score** (first derivative):

$$\frac{\partial}{\partial\mu}\log f = \frac{x_i-\mu}{\sigma^2}$$

**Second derivative** is a constant:

$$\frac{\partial^2}{\partial\mu^2}\log f = -\frac{1}{\sigma^2}$$

**Fisher information** — negate and take expectation, then multiply by $n$:

$$I_1(\mu) = \frac{1}{\sigma^2} \quad\Longrightarrow\quad \boxed{I_n(\mu) = \frac{n}{\sigma^2}}$$

**Cramér-Rao bound:**

$$\text{Var}(\hat\mu) \ge \frac{1}{I_n(\mu)} = \frac{\sigma^2}{n}$$

The sample mean $\bar X$ has variance **exactly** $\sigma^2/n$, so it attains the bound — $\bar X$ is the **efficient** estimator of $\mu$. This is why we just use the sample mean for average return; no unbiased estimator does better.

## Variant — estimating the variance

Now take $\theta = \sigma^2$ (with $\mu$ known). The same procedure gives:

$$I_n(\sigma^2) = \frac{n}{2\sigma^4} \quad\Longrightarrow\quad \text{Var}(\widehat{\sigma^2}) \ge \frac{2\sigma^4}{n}$$

Note the denominator is $\sigma^4$, not $\sigma^2$: estimating variance is "harder" than estimating the mean, and precision degrades faster as the true volatility rises. This is why volatility estimates are relatively stable in calm regimes and noisy in turbulent ones.

## Quant connection

The bound $\sigma^2/n$ is exactly the **squared standard error** that recurs in [[hypothesis_testing]]. Tie it to the Sharpe ratio: testing whether mean return $> 0$ has signal-to-noise

$$\frac{\mu}{\sigma/\sqrt{n}} = \underbrace{\frac{\mu}{\sigma}}_{\text{per-period Sharpe}}\cdot\sqrt{n}$$

More Fisher information (smaller $\sigma$, larger $n$) → more precise $\hat\mu$ → easier to separate a real alpha from noise. **Information, estimation precision, statistical power, and Sharpe are different faces of the same quantity.**

## Related Notes

- [[hypothesis_testing]]
- [[distribution]]
