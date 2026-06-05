#tags: basic-math, problem-set, portfolio, sharpe-ratio, optimization, linear-algebra

# Problem Set — Basic Maths

## Contents

- [[#Problem 1: What is the best weight to combine multiple portfolios?]]

---

## Problem 1: What is the best weight to combine multiple portfolios?

**Problem:** You have $n$ portfolios (or strategies/signals), each with known excess return and standard deviation, and known pairwise correlations. What weights $\mathbf{w}$ combine them into a single portfolio with the **highest Sharpe ratio**?

**Setup:** Let $\boldsymbol{\tilde\mu} = (\tilde\mu_1, \dots, \tilde\mu_n)^\top$ be the vector of **excess** returns ($\tilde\mu_i = \mu_i - r_f$) and $\Sigma$ the $n \times n$ covariance matrix of the portfolios. A combined portfolio with weights $\mathbf{w}$ has excess return $\mathbf{w}^\top\boldsymbol{\tilde\mu}$ and variance $\mathbf{w}^\top\Sigma\mathbf{w}$ (see [[variance_covariance_correlation]]), so its Sharpe ratio is:

$$\text{SR}(\mathbf{w}) = \frac{\mathbf{w}^\top\boldsymbol{\tilde\mu}}{\sqrt{\mathbf{w}^\top\Sigma\mathbf{w}}}.$$

**Solution:**

**Step 1 — Exploit scale-invariance.** SR is unchanged under $\mathbf{w} \to c\mathbf{w}$ (numerator and denominator both scale by $c$). So the optimum is a *direction*, not a magnitude — we are free to impose one normalization. Fix the excess return at 1:

$$\mathbf{w}^\top\boldsymbol{\tilde\mu} = 1.$$

**Step 2 — Reframe as minimum-variance.** With the numerator fixed at 1, maximizing the ratio is equivalent to minimizing the denominator. The problem becomes a convex quadratic program:

$$\min_{\mathbf{w}}\ \mathbf{w}^\top\Sigma\mathbf{w} \quad\text{s.t.}\quad \mathbf{w}^\top\boldsymbol{\tilde\mu} = 1.$$

(Convex because $\Sigma \succeq 0$, so any stationary point is the global optimum.)

**Step 3 — Lagrange.** Form the Lagrangian:

$$\mathcal{L}(\mathbf{w}, \lambda) = \mathbf{w}^\top\Sigma\mathbf{w} - \lambda\big(\mathbf{w}^\top\boldsymbol{\tilde\mu} - 1\big).$$

Stationarity in $\mathbf{w}$ (using $\tfrac{\partial}{\partial\mathbf{w}}\mathbf{w}^\top\Sigma\mathbf{w} = 2\Sigma\mathbf{w}$):

$$2\Sigma\mathbf{w} - \lambda\boldsymbol{\tilde\mu} = 0 \quad\Longrightarrow\quad \boxed{\mathbf{w} \propto \Sigma^{-1}\boldsymbol{\tilde\mu}}$$

**Answer:** The optimal weights are proportional to $\Sigma^{-1}\boldsymbol{\tilde\mu}$ — the excess-return vector pre-multiplied by the inverse covariance (precision) matrix. This is the **tangency portfolio**.

**Step 4 — Maximized Sharpe ratio.**As SR is scale-invariant, so we are free to assume:

$$\mathbf{w} = \Sigma^{-1}\boldsymbol{\tilde\mu}.$$

Plug straight back into the Sharpe formula. Numerator and denominator both reduce to the same quadratic form (using $\Sigma^{-1\top} = \Sigma^{-1}$ and $\Sigma^{-1}\Sigma\,\Sigma^{-1} = \Sigma^{-1}$):

$$\mathbf{w}^\top\boldsymbol{\tilde\mu} = \boldsymbol{\tilde\mu}^\top\Sigma^{-1}\boldsymbol{\tilde\mu}, \qquad \mathbf{w}^\top\Sigma\mathbf{w} = \boldsymbol{\tilde\mu}^\top\Sigma^{-1}\Sigma\,\Sigma^{-1}\boldsymbol{\tilde\mu} = \boldsymbol{\tilde\mu}^\top\Sigma^{-1}\boldsymbol{\tilde\mu}.$$

Therefore

$$\text{SR}_{\max} = \frac{\boldsymbol{\tilde\mu}^\top\Sigma^{-1}\boldsymbol{\tilde\mu}}{\sqrt{\boldsymbol{\tilde\mu}^\top\Sigma^{-1}\boldsymbol{\tilde\mu}}} = \boxed{\sqrt{\boldsymbol{\tilde\mu}^\top\Sigma^{-1}\boldsymbol{\tilde\mu}}}.$$

If there's 

**Step 5 — Decompose into Sharpes and correlations.** Separate the covariance into vols and correlations: $\Sigma = DCD$ with $D = \text{diag}(\sigma_1,\dots,\sigma_n)$ and $C$ the correlation matrix. Write the excess returns as $\boldsymbol{\tilde\mu} = D\mathbf{s}$, where $\mathbf{s} = (\tilde\mu_1/\sigma_1,\dots,\tilde\mu_n/\sigma_n)^\top$ is the vector of **individual Sharpe ratios**. Substituting, the diagonal $D$ factors cancel against $\Sigma^{-1} = D^{-1}C^{-1}D^{-1}$:

$$\text{SR}^2_{\max} = \boldsymbol{\tilde\mu}^\top\Sigma^{-1}\boldsymbol{\tilde\mu} = (D\mathbf{s})^\top D^{-1}C^{-1}D^{-1}(D\mathbf{s}) = \boxed{\mathbf{s}^\top C^{-1}\mathbf{s}}.$$

So the maximized squared Sharpe depends **only on the individual Sharpes and their correlation matrix, the level of return or volatility does not matter

For $n=2$, with individual Sharpes $s_1, s_2$ and correlation $\rho$, $C^{-1} = \frac{1}{1-\rho^2}\begin{psmallmatrix} 1 & -\rho \\ -\rho & 1\end{psmallmatrix}$ gives the explicit closed form:

$$\text{SR}^2_{\max} = \frac{s_1^2 - 2\rho\, s_1 s_2 + s_2^2}{1-\rho^2}.$$


**Special case — equal individual Sharpes.** If every individual Sharpe is 1 ($\mathbf{s} = \mathbf 1$), then $\text{SR}^2_{\max} = \mathbf{1}^\top C^{-1}\mathbf{1}$ (the sum of all entries of $C^{-1}$):

- **Uncorrelated** ($C = I$): $\text{SR}_{\max} = \sqrt{n}$ — the $\sqrt{n}$ diversification law; independent Sharpes add in quadrature.
- **Bounds** ($n=2$, correlation $\rho$): $\text{SR}^2_{\max} = \frac{2}{1+\rho}$. As $\rho\to+1$ (redundant), $\text{SR}_{\max}\to 1$ — never below the best single stream. As $\rho\to-1$ (perfect hedge), $\text{SR}_{\max}\to\infty$ — the long–short combination cancels risk while returns add. So $\text{SR}_{\max}\in[1,\infty)$.
- **Caveat:** the unbounded upside is an artifact of known $\boldsymbol{\tilde\mu},\Sigma$. Near $\rho=\pm1$, $C$ is near-singular and $\Sigma^{-1}$ amplifies estimation error along exactly those directions — realized Sharpe collapses. See the estimation-error caveat below.

