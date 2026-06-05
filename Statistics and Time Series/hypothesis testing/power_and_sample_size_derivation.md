#tags: hypothesis-testing, power, type-I-error, type-II-error, sample-size, z-test

# Power, Type I/II Error

## Goal

Derive the relationship linking $\alpha$, $\beta$, effect size, and sample size:

$$z_{1-\alpha} + z_{1-\beta} = \frac{\mu_1 - \mu_0}{\sigma/\sqrt{n}}$$

starting from the sampling distribution of the mean. No formula here is memorized — each falls out of the previous step.

## Setup

One-sided z-test, known variance, simple (point) alternative:

$$H_0: \mu = \mu_0 \qquad H_1: \mu = \mu_1 \quad (\mu_1 > \mu_0)$$

Sample size $n$, known population SD $\sigma$. Reject $H_0$ when $\bar{X} > c$ for some critical value $c$.

## Foundational fact

The sampling distribution of the mean, and its standardization:

$$\bar{X} \sim \mathcal{N}\!\left(\mu, \frac{\sigma^2}{n}\right) \quad\Longrightarrow\quad Z = \frac{\bar{X} - \mu}{\sigma/\sqrt{n}} \sim \mathcal{N}(0,1)$$

The quantity $\sigma/\sqrt{n}$ is the **standard error (SE)**. Its $\sqrt{n}$ scaling is the source of every sample-size effect downstream.

## Step 1 — Critical value from $\alpha$ (work under $H_0$)

$\alpha$ is the false-positive rate, computed **assuming the null is true**:

$$\alpha = P(\text{reject} \mid H_0) = P(\bar{X} > c \mid \mu = \mu_0)$$

Standardize using $\mu_0$ (we are in the null world):

$$\alpha = P\!\left(Z > \frac{c-\mu_0}{\sigma/\sqrt{n}}\right)$$

For this tail area to equal $\alpha$, the cutoff in $Z$-space must be the quantile $z_{1-\alpha}$ (**quantile inversion**):

$$\frac{c-\mu_0}{\sigma/\sqrt{n}} = z_{1-\alpha} \quad\Longrightarrow\quad \boxed{c = \mu_0 + z_{1-\alpha}\frac{\sigma}{\sqrt{n}}}$$

## Step 2 — $\beta$ from the critical value (work under $H_1$)

$\beta$ is computed in a **different world** — the alternative. Same cutoff $c$, but the true mean is now $\mu_1$:

$$\beta = P(\text{fail to reject} \mid H_1) = P(\bar{X} < c \mid \mu = \mu_1)$$

Standardize again, this time subtracting $\mu_1$ (the true mean now):

$$\beta = \Phi\!\left(\frac{c-\mu_1}{\sigma/\sqrt{n}}\right)$$

## Step 3 — Substitute $c$ and simplify

Plug the $c$ from Step 1 into the $\beta$ expression and split the fraction:

$$\frac{c-\mu_1}{\sigma/\sqrt{n}} = \frac{\mu_0 + z_{1-\alpha}\frac{\sigma}{\sqrt{n}} - \mu_1}{\sigma/\sqrt{n}} = z_{1-\alpha} + \frac{\mu_0-\mu_1}{\sigma/\sqrt{n}} = z_{1-\alpha} - \frac{\mu_1-\mu_0}{\sigma/\sqrt{n}}$$

So:

$$\boxed{\beta = \Phi\!\left(z_{1-\alpha} - \frac{\mu_1-\mu_0}{\sigma/\sqrt{n}}\right)}$$

The term $\dfrac{\mu_1-\mu_0}{\sigma/\sqrt{n}}$ is the **standardized effect size** — the distance between the two hypotheses measured in standard errors.

## Step 4 — Symmetric form via one more inversion

Apply $\Phi^{-1}$ to both sides. Use the **normal symmetry identity** $\Phi^{-1}(\beta) = -z_{1-\beta}$ (i.e. $z_p = -z_{1-p}$ — this is the most common place to slip a sign):

$$-z_{1-\beta} = z_{1-\alpha} - \frac{\mu_1-\mu_0}{\sigma/\sqrt{n}}$$

Rearrange:

$$\boxed{z_{1-\alpha} + z_{1-\beta} = \frac{\mu_1 - \mu_0}{\sigma/\sqrt{n}}}$$

**Punchline:** the sum of the two critical z-values equals the standardized effect size.

## Sample-size formula (corollary)

Solve the symmetric form for $n$:

$$n = \frac{(z_{1-\alpha} + z_{1-\beta})^2 \,\sigma^2}{(\mu_1 - \mu_0)^2}$$

## The reusable toolkit

Every step above is one of five recurring moves:

1. **Standardization** $Z = (\bar{X}-\mu)/(\sigma/\sqrt{n})$ — used twice (once under $H_0$, once under $H_1$).
2. **Quantile inversion** — "tail area $= p$" $\Leftrightarrow$ "cutoff $= z_p$", via $\Phi$ / $\Phi^{-1}$.
3. **Normal symmetry** $z_p = -z_{1-p}$ — collapses the messy form into the clean symmetric one.
4. **Standard error $\sigma/\sqrt{n}$** and its $\sqrt{n}$ scaling — source of all sample-size effects.
5. **Switching worlds** — $\alpha$ is computed under $H_0$, $\beta$ under $H_1$. Same data, same cutoff, two assumed truths. Keep straight which $\mu$ you subtract at each step.

## Caveats

- Requires a **simple (point) alternative**. For composite $H_1: \mu > \mu_0$, $\beta$ becomes a *function* of the true $\mu$ (the power curve), not a single number.
- Exact only for **known $\sigma$ and normality**. For the t-test (unknown $\sigma$), the same skeleton holds but the alternative distribution is *noncentral* t with noncentrality parameter $\delta\sqrt{n}/\sigma$ — no equally clean closed form.

## Related Notes

- [[type_I_and_II_errors]]
- [[p_value]]
