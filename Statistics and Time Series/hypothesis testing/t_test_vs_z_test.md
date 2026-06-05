#tags: hypothesis-testing, z-test, t-test, t-distribution, sample-size, noncentral-t

# z-test vs. t-test

## The one-line distinction

Use a **z-test** when the population SD $\sigma$ is **known**; use a **t-test** when $\sigma$ is **unknown** and estimated from the sample with $s$. In practice $\sigma$ is almost never known, so the t-test is the real-world default.

## The statistics

$$z = \frac{\bar{X} - \mu_0}{\sigma/\sqrt{n}} \sim \mathcal{N}(0,1) \qquad t = \frac{\bar{X} - \mu_0}{s/\sqrt{n}} \sim t_{\,n-1}$$

Same shape, one difference: $\sigma \to s$ in the denominator.

## Why estimating $\sigma$ changes the distribution

$s$ is itself a **random variable** — it wobbles from sample to sample. So the t-statistic has randomness in **both** numerator and denominator, whereas z has it only in the numerator (fixed $\sigma$ below).

That extra denominator noise makes large values more likely: sometimes $s$ underestimates $\sigma$ by chance, inflating the ratio. Hence the t-distribution has **fatter tails** than the normal, so its critical values are pushed further out — the t-test is **more conservative**. That is the price of not knowing $\sigma$.

| degree of freedom ($\nu$) | two-sided 5% critical value |
| ------------------------- | --------------------------- |
| 10                        | 2.23                        |
| 25                        | 2.06                        |
| 100                       | 1.98                        |
| $\infty$ (z)              | 1.96                        |

in one-sample test, degree of freedom df = n - 1.
## Effect on power / Type II error

For the **null side** (critical value, Type I error), moving from z to t is a clean swap: $\sigma \to s$, $z_{1-\alpha} \to t_{1-\alpha,\,n-1}$.

The **power side does not swap cleanly.** Under $H_0$ the statistic is **central** $t_{n-1}$; under $H_1$ it is **noncentral** $t_{n-1}(\delta)$ with noncentrality $\delta = (\mu_1-\mu_0)/(\sigma/\sqrt{n})$. These are *different shapes*, not one shifted, so the symmetric identity $z_{1-\alpha}+z_{1-\beta}=\delta$ has **no t-analogue** — power and sample size are computed **numerically**.

| | Under $H_0$ | Under $H_1$ | Closed-form power? |
|---|---|---|---|
| z-test | $N(0,1)$ | $N(\delta, 1)$ — same shape, shifted | Yes |
| t-test | central $t_\nu$ | noncentral $t_\nu(\delta)$ — different shape | No |

See [[power_and_sample_size_derivation]] for the clean z-test derivation this contrasts with.

## Quick decision guide

- $\sigma$ known → z-test (rare in practice).
- $\sigma$ unknown, $n$ large → t-test 

## Related Notes

- [[power_and_sample_size_derivation]]
