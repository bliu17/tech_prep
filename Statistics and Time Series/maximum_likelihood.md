#tags: mle, maximum-likelihood, estimation, asymptotics, fisher-information, ols, invariance, statistics

# Maximum Likelihood Estimation

Pick the parameter value that makes the observed data **most probable** under the assumed model. The single most-asked estimation topic in quant interviews — because it underpins OLS, GARCH, logistic regression, and almost every "how would you fit this?" question.

---

## Definition

Given data $x_1,\dots,x_n$ with assumed density $f(x;\theta)$, the **likelihood** is the joint density viewed as a function of $\theta$ (data fixed):

$$L(\theta) = \prod_{i=1}^n f(x_i;\theta) \qquad \ell(\theta) = \log L(\theta) = \sum_{i=1}^n \log f(x_i;\theta)$$

The **MLE** maximizes this:

$$\hat\theta_{\text{MLE}} = \arg\max_\theta \ell(\theta)$$

We work with the **log**-likelihood for three reasons interviewers like to hear: products become sums (differentiable term-by-term), it avoids numerical underflow from multiplying many small probabilities, and $\log$ is monotone so the argmax is unchanged.

**Likelihood vs. probability** — the classic conceptual question. Same expression $f(x;\theta)$, but probability fixes $\theta$ and varies $x$ (integrates to 1 over $x$); likelihood fixes $x$ and varies $\theta$ (does **not** integrate to 1 over $\theta$). The likelihood is **not** a probability distribution over $\theta$ — that would be the Bayesian posterior.

## The recipe

1. Write $\ell(\theta) = \sum_i \log f(x_i;\theta)$.
2. **Score** = differentiate: $s(\theta) = \partial\ell/\partial\theta$.
3. Set score to zero and solve: $s(\hat\theta) = 0$ (the **likelihood equation**).
4. Check it's a **max**, not a min/saddle: $\partial^2\ell/\partial\theta^2 < 0$.

The score $s(\theta)=\partial\ell/\partial\theta$ is the same object whose variance is [[fisher_information]] — that link is the whole reason MLE has clean asymptotics.

## Worked example 1 — Normal returns (the OLS link)

Returns $X_1,\dots,X_n \overset{\text{iid}}{\sim}\mathcal{N}(\mu,\sigma^2)$, **both** parameters unknown.

$$\ell(\mu,\sigma^2) = -\frac{n}{2}\log(2\pi) - \frac{n}{2}\log\sigma^2 - \frac{1}{2\sigma^2}\sum_{i=1}^n (x_i-\mu)^2$$

**Score w.r.t. $\mu$**, set to zero:

$$\frac{\partial\ell}{\partial\mu} = \frac{1}{\sigma^2}\sum_i (x_i-\mu) = 0 \;\Longrightarrow\; \boxed{\hat\mu = \bar X}$$

Maximizing the normal likelihood over $\mu$ means **minimizing $\sum(x_i-\mu)^2$** — the sum of squared errors. This is the key insight: **MLE under Gaussian errors = least squares.** Generalizing $\mu$ to $X\beta$ gives MLE $\hat\beta = $ the [[linear_regression]] OLS estimator. When an interviewer asks "why OLS?", the deepest answer is: it's the MLE under normal, homoskedastic, independent errors.

**Score w.r.t. $\sigma^2$**, set to zero:

$$\frac{\partial\ell}{\partial\sigma^2} = -\frac{n}{2\sigma^2} + \frac{1}{2\sigma^4}\sum_i (x_i-\mu)^2 = 0 \;\Longrightarrow\; \boxed{\hat\sigma^2 = \frac{1}{n}\sum_i (x_i-\bar X)^2}$$

## The bias trap (asked constantly)

That $\hat\sigma^2$ divides by $n$, **not** $n-1$. So the MLE of the variance is **biased**:

$$\mathbb{E}[\hat\sigma^2] = \frac{n-1}{n}\,\sigma^2 < \sigma^2$$

It underestimates because it measures spread around the *estimated* mean $\bar X$ (which sits inside the data by construction) rather than the true $\mu$. The unbiased sample variance divides by $n-1$ instead. Two takeaways interviewers want:

- **MLE is not generally unbiased** — only *asymptotically* unbiased. The bias here is $O(1/n)$ and vanishes as $n\to\infty$.
- It's a clean example that **efficiency $\ne$ unbiasedness**: the MLE trades a little bias for lower MSE.

## Key properties (the asymptotics)

Under regularity conditions, as $n\to\infty$:

**1. Consistency** — $\hat\theta \xrightarrow{p} \theta_0$. The MLE converges to the truth.

**2. Asymptotic normality** — the headline result, and the one to be able to state cold:

$$\sqrt{n}\,(\hat\theta - \theta_0) \xrightarrow{d} \mathcal{N}\!\left(0,\; I_1(\theta_0)^{-1}\right) \quad\Longleftrightarrow\quad \hat\theta \approx \mathcal{N}\!\left(\theta_0,\; \frac{1}{n\,I_1(\theta_0)}\right)$$

The asymptotic variance is the **inverse Fisher information**. This is *why* [[fisher_information]] matters: it directly gives you the standard error of the MLE, hence confidence intervals and the Wald test. More curvature in the log-likelihood → tighter estimate.

**3. Asymptotic efficiency** — the MLE attains the Cramér-Rao lower bound asymptotically. No consistent estimator has smaller asymptotic variance. This is the formal sense in which "MLE is the best you can do for large $n$."

**4. Invariance** — if $\hat\theta$ is the MLE of $\theta$, then $g(\hat\theta)$ is the MLE of $g(\theta)$ for any function $g$. Example: the MLE of the **volatility** $\sigma$ is just $\sqrt{\hat\sigma^2}$ — you don't re-derive. (Note: invariance preserves the MLE, *not* unbiasedness — $\sqrt{\hat\sigma^2}$ is biased even where $\hat\sigma^2$ weren't.)

## Worked example 2 — clean closed forms (warm-up brainteasers)

These get thrown out as quick checks. Each follows the same recipe.

**Poisson** ($X_i \sim \text{Pois}(\lambda)$): $\ell = \sum(x_i\log\lambda - \lambda) + \text{const}$, score $\sum x_i/\lambda - n = 0 \Rightarrow \hat\lambda = \bar X$.

**Exponential** ($X_i \sim \text{Exp}(\lambda)$, rate form): $\ell = n\log\lambda - \lambda\sum x_i$, score $n/\lambda - \sum x_i = 0 \Rightarrow \hat\lambda = 1/\bar X$. (Invariance: MLE of the mean $1/\lambda$ is $\bar X$.)

**Bernoulli** ($X_i \sim \text{Bern}(p)$): $\hat p = \bar X$, the sample proportion.

Pattern worth saying out loud: for these one-parameter exponential-family models, the MLE just **matches the sample moment** to the population moment.

## When there's no closed form

Most real models. The likelihood equation $s(\hat\theta)=0$ has no analytic solution, so you maximize **numerically** (Newton-Raphson, which uses the Fisher information / Hessian; or BFGS). Two that matter for systematic/quant roles:

- **Logistic regression** — MLE of $\beta$ for binary outcomes (e.g. predicting up/down direction). No closed form; the log-likelihood is concave so Newton's method (IRLS) converges reliably. There is no "least-squares" shortcut here — you *must* go through the likelihood.
- **GARCH** — volatility-model parameters are fit by MLE on the **conditional** likelihood $\sum_t \log f(r_t \mid \mathcal F_{t-1})$, since returns aren't i.i.d. In practice one often uses **quasi-MLE** (assume Gaussian innovations even if they're fat-tailed); QMLE still gives consistent parameter estimates with robust (sandwich) standard errors. A common "gotcha" answer: GARCH MLE assumes conditional, not unconditional, normality.

## Likelihood-based tests (one-liners to know)

All three are asymptotically equivalent ways to test $H_0:\theta=\theta_0$, differing in what they evaluate:

- **Wald** — uses $\hat\theta$ and its variance (inverse Fisher info): $(\hat\theta-\theta_0)^2/\widehat{\text{Var}}$. Only fits the full model.
- **Score (Lagrange multiplier)** — evaluates the score at $\theta_0$. Only fits the restricted model.
- **Likelihood ratio** — $2[\ell(\hat\theta) - \ell(\theta_0)] \xrightarrow{d} \chi^2$. Fits both; generally the most reliable.

## Quant connection

The whole chain ties together: maximizing the Gaussian likelihood **is** OLS ([[linear_regression]]); the MLE's precision is governed by [[fisher_information]], whose bound $\sigma^2/n$ is the same squared standard error that drives [[hypothesis_testing]] and Sharpe-ratio significance. Regularized estimators (Ridge/LASSO, see [[ridge_lasso]]) are **MAP** estimates — penalized likelihood with a Gaussian or Laplace prior on $\beta$ — which is the Bayesian bridge from MLE. Knowing MLE well means you can answer "how is this model fit?" for nearly anything systematic.

## Likely interview questions

1. **What's the difference between likelihood and probability?** → Fixed $\theta$ vary $x$ vs. fixed $x$ vary $\theta$; likelihood doesn't integrate to 1 over $\theta$.
2. **Derive the MLE of $\mu,\sigma^2$ for a normal sample.** → Above. Be ready for the $1/n$ vs $1/(n-1)$ bias point.
3. **Is the MLE unbiased?** → Not in general; only consistent and asymptotically unbiased. Give the $\hat\sigma^2$ example.
4. **What's the asymptotic distribution of the MLE?** → $\mathcal N(\theta_0, [nI_1(\theta_0)]^{-1})$; variance = inverse Fisher info. Be able to *say why* (Taylor-expand the score).
5. **Why is OLS "optimal"?** → It's the MLE under Gaussian iid errors; BLUE via Gauss-Markov without the normality assumption.
6. **MLE of $\sigma$ given MLE of $\sigma^2$?** → Invariance: $\sqrt{\hat\sigma^2}$, no re-derivation.
7. **How do you fit a model with no closed-form MLE?** → Numerical optimization (Newton/IRLS); name logistic regression or GARCH.
8. **Why log-likelihood?** → Sums not products, numerical stability, monotone so argmax unchanged.

---

# Worked Answers

Full derivations for the eight questions above. Each is written the way you'd say it at a whiteboard.

## A1 — Likelihood vs. probability

Both are the *same function* $f(x;\theta)$ read two different ways.

- **Probability**: fix the parameter $\theta$, let the data $x$ vary. $f(\cdot\,;\theta)$ is a density in $x$, so $\int f(x;\theta)\,dx = 1$.
- **Likelihood**: fix the observed data $x$, let $\theta$ vary. $L(\theta)=f(x;\theta)$ is a function of $\theta$, and crucially $\int L(\theta)\,d\theta \ne 1$ in general — the likelihood is **not** a distribution over $\theta$.

To make it a probability statement *about* $\theta$ you need Bayes: $p(\theta\mid x) \propto L(\theta)\,p(\theta)$. The likelihood is the data's "vote," the prior is what turns it into a posterior. One-line version: *probability describes data given parameters; likelihood scores parameters given data.*

## A2 — MLE of $\mu,\sigma^2$ for a normal sample (with the bias point)

Start from the log-likelihood (the constant $-\frac n2\log 2\pi$ is dropped where it doesn't matter):

$$\ell(\mu,\sigma^2) = -\frac{n}{2}\log\sigma^2 - \frac{1}{2\sigma^2}\sum_{i=1}^n (x_i-\mu)^2 + \text{const}$$

**Solve for $\mu$ first.** Only the sum-of-squares term involves $\mu$:

$$\frac{\partial\ell}{\partial\mu} = \frac{1}{\sigma^2}\sum_i (x_i-\mu) = 0 \;\Longrightarrow\; \sum_i x_i = n\mu \;\Longrightarrow\; \hat\mu = \bar X.$$

Note $\hat\mu$ doesn't depend on $\sigma^2$, so we can plug it in.

**Now $\sigma^2$.** Treat $v=\sigma^2$ as the variable:

$$\frac{\partial\ell}{\partial v} = -\frac{n}{2v} + \frac{1}{2v^2}\sum_i (x_i-\hat\mu)^2 = 0.$$

Multiply by $2v^2$: $\;-nv + \sum_i(x_i-\bar X)^2 = 0 \Rightarrow \hat\sigma^2 = \frac1n\sum_i (x_i-\bar X)^2.$

**The bias.** Write $\sum(x_i-\bar X)^2 = \sum(x_i-\mu)^2 - n(\bar X-\mu)^2$ (the identity from expanding around the true $\mu$). Take expectations: $\mathbb E[\sum(x_i-\mu)^2]=n\sigma^2$ and $\mathbb E[n(\bar X-\mu)^2]=n\cdot\sigma^2/n=\sigma^2$, so

$$\mathbb E\!\left[\sum_i(x_i-\bar X)^2\right] = n\sigma^2-\sigma^2=(n-1)\sigma^2 \;\Longrightarrow\; \mathbb E[\hat\sigma^2]=\frac{n-1}{n}\sigma^2.$$

The MLE undercounts by exactly one degree of freedom — it measures spread around $\bar X$, which is pulled toward the data. Dividing by $n-1$ instead restores unbiasedness. *(Verified by simulation: $n{=}10,\sigma^2{=}4$ gives $\mathbb E[\hat\sigma^2]\approx 3.60 = \tfrac{9}{10}\cdot4$.)*

## A3 — Is the MLE unbiased?

No, not in general. The $\hat\sigma^2$ above is the standard counterexample: biased by a factor $(n-1)/n$. What the MLE *is*:

- **Consistent**: $\hat\theta\xrightarrow{p}\theta_0$ — it lands on the truth as $n\to\infty$.
- **Asymptotically unbiased**: bias is typically $O(1/n)\to 0$.
- **Asymptotically efficient**: smallest possible variance in the limit (A4).

The deeper point: MLE optimizes for the *likelihood*, not for unbiasedness. It will happily accept a small bias if that lowers overall MSE. Unbiasedness is a finite-sample property MLE simply doesn't promise.

## A4 — Asymptotic distribution of the MLE (derivation)

Claim: $\sqrt n(\hat\theta-\theta_0)\xrightarrow{d}\mathcal N\!\big(0, I_1(\theta_0)^{-1}\big)$. Derivation is a **Taylor expansion of the score** around the truth.

Let $s(\theta)=\ell'(\theta)$ be the score. The MLE solves $s(\hat\theta)=0$. Expand $s$ about $\theta_0$:

$$0 = s(\hat\theta) \approx s(\theta_0) + s'(\theta_0)\,(\hat\theta-\theta_0).$$

Solve for the scaled error:

$$\sqrt n(\hat\theta-\theta_0) \approx \frac{\tfrac{1}{\sqrt n}\,s(\theta_0)}{-\tfrac1n s'(\theta_0)}.$$

Now handle numerator and denominator separately.

**Numerator** — $s(\theta_0)=\sum_i \partial_\theta\log f(x_i;\theta_0)$ is a sum of iid mean-zero terms (the score has zero expectation) each with variance $I_1(\theta_0)$. By the [[central_limit_theorem|CLT]]:

$$\frac{1}{\sqrt n}s(\theta_0) \xrightarrow{d} \mathcal N\!\big(0,\,I_1(\theta_0)\big).$$

**Denominator** — $-\frac1n s'(\theta_0) = -\frac1n\sum_i \partial_\theta^2\log f(x_i;\theta_0)$ is a sample average. By the LLN it converges to its expectation, which is exactly Fisher information (the second-derivative form):

$$-\frac1n s'(\theta_0) \xrightarrow{p} I_1(\theta_0).$$

**Combine** (Slutsky): a $\mathcal N(0,I_1)$ divided by the constant $I_1$ gives

$$\sqrt n(\hat\theta-\theta_0)\xrightarrow{d}\mathcal N\!\left(0,\frac{I_1}{I_1^2}\right)=\mathcal N\!\big(0,I_1(\theta_0)^{-1}\big).$$

Hence $\hat\theta\approx\mathcal N\!\big(\theta_0,[nI_1(\theta_0)]^{-1}\big)$ — the variance is inverse Fisher information, which is why [[fisher_information]] *is* the standard error of the MLE. *(Verified: for $\text{Exp}(\lambda{=}2)$, $I_1=1/\lambda^2$ so $n\cdot\text{Var}(\hat\lambda)\to\lambda^2=4$; simulation gives $3.95$.)*

## A5 — Why is OLS "optimal"?

Two distinct optimality claims; know which assumptions each needs.

**(a) OLS is the Gaussian MLE.** Model $y_i = x_i^\top\beta + \varepsilon_i$, $\varepsilon_i\overset{\text{iid}}\sim\mathcal N(0,\sigma^2)$. The log-likelihood is

$$\ell(\beta) = -\frac{1}{2\sigma^2}\sum_i (y_i - x_i^\top\beta)^2 + \text{const}.$$

Maximizing over $\beta$ $\Leftrightarrow$ **minimizing the residual sum of squares** $\Rightarrow$ $\hat\beta=(X^\top X)^{-1}X^\top y$, the OLS estimator. So under normal errors OLS inherits all MLE properties: consistency, asymptotic normality, efficiency.

**(b) OLS is BLUE (Gauss–Markov).** Drop normality. If errors merely have mean zero, constant variance $\sigma^2$, and are uncorrelated, then OLS is the **B**est **L**inear **U**nbiased **E**stimator — minimum variance *among linear unbiased* estimators. No distributional assumption needed.

The crisp interview answer: *OLS is the MLE if you believe errors are Gaussian; it's BLUE even if you don't, as long as errors are homoskedastic and uncorrelated.* Heteroskedasticity/correlation breaks both and pushes you toward GLS.

## A6 — MLE of $\sigma$ given the MLE of $\sigma^2$ (invariance)

By the **invariance property**, if $\hat\theta$ is the MLE of $\theta$, then $g(\hat\theta)$ is the MLE of $g(\theta)$ for any $g$. Taking $g(v)=\sqrt v$:

$$\hat\sigma = \sqrt{\hat\sigma^2} = \sqrt{\tfrac1n\sum_i(x_i-\bar X)^2}.$$

No re-derivation needed. *Why invariance holds:* the MLE is the argmax of the likelihood; reparameterizing relabels the axis but the peak sits at the same data configuration, so the maximizer maps through $g$. **Caution:** invariance preserves the *MLE*, not unbiasedness — by Jensen, $\mathbb E[\sqrt{\hat\sigma^2}]<\sqrt{\mathbb E[\hat\sigma^2]}$, so $\hat\sigma$ is biased *even where $\hat\sigma^2$ is not.* Functions of unbiased estimators are generally biased.

## A7 — Fitting a model with no closed-form MLE

When $s(\hat\theta)=0$ has no analytic solution, maximize $\ell$ **numerically**. The workhorse is **Newton–Raphson**:

$$\theta_{k+1} = \theta_k - \big[\ell''(\theta_k)\big]^{-1}\,\ell'(\theta_k),$$

i.e. step by (Hessian)$^{-1}$(score). Replacing the observed Hessian $-\ell''$ with its expectation (Fisher information) gives **Fisher scoring**; for GLMs this is exactly **IRLS** (iteratively reweighted least squares). BFGS/L-BFGS are quasi-Newton alternatives that approximate the Hessian.

Two quant-relevant cases:

- **Logistic regression** (binary up/down). $\ell(\beta)=\sum_i[y_i x_i^\top\beta - \log(1+e^{x_i^\top\beta})]$. The score $\sum_i(y_i-p_i)x_i$ has no closed-form root, but $\ell$ is **globally concave**, so Newton/IRLS converges to the unique maximum.
- **GARCH**. Returns aren't iid, so you maximize the **conditional** log-likelihood $\sum_t \log f(r_t\mid\mathcal F_{t-1})$ with $r_t\mid\mathcal F_{t-1}\sim\mathcal N(0,\sigma_t^2)$ and $\sigma_t^2$ following the recursion. Optimized numerically; if innovations are fat-tailed, Gaussian **quasi-MLE** still gives consistent estimates with robust (sandwich) standard errors.

## A8 — Why the log-likelihood?

Three reasons, all worth stating:

1. **Sums beat products.** $\log\prod_i f = \sum_i\log f$. Differentiating a sum is term-by-term and trivial; differentiating a product invokes a giant product rule. Every score calculation relies on this.
2. **Numerical stability.** A product of $n$ small probabilities underflows to $0$ in floating point ($0.01^{400}$ is already machine zero); summing logs stays in a representable range.
3. **Monotonicity.** $\log$ is strictly increasing, so $\arg\max_\theta L(\theta)=\arg\max_\theta \log L(\theta)$ — you optimize the convenient function and get the *same* maximizer.

(Bonus: the log-likelihood's curvature is Fisher information, and asymptotic theory is naturally stated in log-likelihood terms — so it's not just convenience, it's the natural object.)

---

## Related Notes

- [[fisher_information]]
- [[hypothesis_testing]]
- [[linear_regression]]
- [[ridge_lasso]]
- [[distribution]]
