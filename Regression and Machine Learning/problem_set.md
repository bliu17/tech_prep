#tags: regression, OLS, problem-set, multicollinearity, linear-algebra

# Problem Set — Regression

## Contents

- [[#Problem 1: Assume you regress y ~ x1 + x2, how does changing correlation of x1 and x2 impact the Variance of OLS Estimator]]
- [[#Problem 2: What is the geometric meaning of OLS?]]
- [[#Problem 3: How many points are needed to halve the standard deviation of beta?]]
- [[#Problem 4: Is there a closed-form solution to constrained linear regression?]]
- [[#Problem 5: If you duplicate the dataset m times, how do beta, t-stats, p-values, and R-squared change?]]

---

## Problem 1: Assume you regress y ~ x1 + x2, how does changing correlation of x1 and x2 impact the Variance of OLS Estimator

**Problem:** Illustrate how the variance of the OLS solution $\hat{\beta}$ changes when the correlation between feature 1 and feature 2 changes, in a two-variable regression.

**Setup:** $n$ observations, 2 features, both demeaned and standardized (unit variance). Let $\rho = \text{Corr}(X_1, X_2)$.

**Solution:**

Since features are standardized, $X^\top X$ has a clean form:

$$X^\top X = n \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix}$$

The diagonal entries are $\sum x_{i1}^2 = \sum x_{i2}^2 = n$ (unit variance), and the off-diagonal is $\sum x_{i1} x_{i2} = n\rho$ (standardized covariance = correlation).

Inverting:

$$(X^\top X)^{-1} = \frac{1}{n(1 - \rho^2)} \begin{bmatrix} 1 & -\rho \\ -\rho & 1 \end{bmatrix}$$

Since the variance covariance matrix of beta $\text{Var}(\hat{\beta}) = \hat{\sigma}^2 (X^\top X)^{-1}$, the variance $\beta$ is:

$$\text{Var}(\hat{\beta}_1) = \text{Var}(\hat{\beta}_2) = \frac{\hat{\sigma}^2}{n(1 - \rho^2)}$$

**Interpretation:**

- As $\rho \to 0$: variance $= \hat{\sigma}^2 / n$ — the baseline, as if the two features were independent.
- As $\rho \to \pm 1$: $(1 - \rho^2) \to 0$, so variance $\to \infty$. The matrix $X^\top X$ becomes nearly singular (determinant $= n^2(1-\rho^2) \to 0$), meaning the columns of $X$ are nearly linearly dependent and OLS cannot separately identify the two coefficients.
- The covariance between $\hat{\beta}_1$ and $\hat{\beta}_2$ is $\frac{-\rho \hat{\sigma}^2}{n(1-\rho^2)}$ — negative when $\rho > 0$, meaning the two estimates are pulled in opposite directions when features are positively correlated.

---

## Problem 2: What is the geometric meaning of OLS when regressing y towards X?

**Problem:** Explain what OLS is minimizing geometrically, and how the solution relates to projection.

**Solution:**
 The solution $\hat{y} = X\hat{\beta}$ is the **orthogonal projection** of $y$ onto the column space of $X$.

For OLS, in terms of geometry, it's minimizing sum of square of vertical distance (not perpendicular distance!) of each point toward the fitted line. If flip the regression, X ~ Y, it becomes minimizing square of horizontal distance of each point to a fitted line.

**Why "ordinary":** OLS minimizes vertical distances in $y$ only , treating $X$ as fixed and noise-free. This contrasts with **total least squares** (orthogonal regression), which minimizes perpendicular (点到直线的垂直距离) distances to the fitted line, and it's assumption is  $X$ itself has measurement error.

---

## Problem 3: How many points are needed to halve the standard deviation of beta of a simple linear regression?

**Problem:** In OLS, the standard deviation of $\hat{\beta}_1$ is $\sigma$. How many data points are needed on average to reduce it to $\sigma/2$?

**Solution:**

Starting from the variance formula:

$$\text{Var}(\hat{\beta}_1) = \frac{\sigma^2}{\sum(x_i - \bar{x})^2}$$

Since $\sum(x_i - \bar{x})^2 = n \cdot \text{SD}(X)^2$:

$$\text{Var}(\hat{\beta}_1) = \frac{\sigma^2}{n \cdot \text{SD}(X)^2} \implies \text{SE}(\hat{\beta}_1) = \frac{\sigma}{\sqrt{n} \cdot \text{SD}(X)}$$

To halve the SE, we need:

$$\frac{\sigma}{\sqrt{n'} \cdot \text{SD}(X)} = \frac{1}{2} \cdot \frac{\sigma}{\sqrt{n} \cdot \text{SD}(X)}$$

$$\sqrt{n'} = 2\sqrt{n} \implies n' = 4n$$

**Answer:** You need **4 times** as many data points. This is the standard $\sqrt{n}$ result — halving the standard error requires quadrupling the sample size.

---

## Problem 4: Is there a closed-form solution to constrained linear regression?

**Problem:** OLS minimizes $\|y - X\beta\|^2$. Now add constraints on $\beta$. Is there a closed-form solution under (a) the equality constraint $A\beta = b$, and (b) the inequality constraints $A\beta \leq b$ or $A\beta \geq b$? Here $A$ is $m \times p$ (with $m < p$, so the constraints don't fully pin down $\beta$) and $b$ is $m \times 1$.

**Setup:** Let $\hat\beta_{OLS} = (X^\top X)^{-1}X^\top y$ be the unconstrained solution. The objective $\|y - X\beta\|^2$ is convex (quadratic, PSD Hessian $2X^\top X$), which is what makes the constrained problems tractable.

### Part (a): Equality constraint $A\beta = b$ — YES, closed form

Minimize $\|y - X\beta\|^2$ subject to $A\beta = b$. Because we have a **convex quadratic objective with linear equality constraints**, Lagrange multipliers give an exact solution. Form the Lagrangian:

$$\mathcal{L}(\beta, \lambda) = (y - X\beta)^\top (y - X\beta) + 2\lambda^\top (A\beta - b)$$

(the factor 2 is cosmetic, to clean up the algebra). The first-order conditions are:

$$\frac{\partial \mathcal{L}}{\partial \beta} = -2X^\top(y - X\beta) + 2A^\top\lambda = 0 \quad\Rightarrow\quad X^\top X\beta = X^\top y - A^\top\lambda$$

$$\frac{\partial \mathcal{L}}{\partial \lambda} = A\beta - b = 0$$

From the first equation, $\beta = \hat\beta_{OLS} - (X^\top X)^{-1}A^\top\lambda$. Substitute into $A\beta = b$ to solve for $\lambda$:

$$A\hat\beta_{OLS} - A(X^\top X)^{-1}A^\top\lambda = b \;\Rightarrow\; \lambda = \big[A(X^\top X)^{-1}A^\top\big]^{-1}\big(A\hat\beta_{OLS} - b\big).$$

Plugging back gives the **closed-form constrained estimator**:

$$\boxed{\hat\beta_C = \hat\beta_{OLS} - (X^\top X)^{-1}A^\top\big[A(X^\top X)^{-1}A^\top\big]^{-1}\big(A\hat\beta_{OLS} - b\big)}$$

**Interpretation:**

- The estimator is the **unconstrained OLS solution, projected onto the constraint surface** $\{ \beta : A\beta = b\}$. The correction term is exactly what's needed to cancel the constraint violation $A\hat\beta_{OLS} - b$.
- If OLS already satisfies the constraint ($A\hat\beta_{OLS} = b$), the correction is zero and $\hat\beta_C = \hat\beta_{OLS}$.
- This is the engine behind testing linear restrictions (the **F-test** compares $\text{SSR}$ of this restricted fit to the unrestricted fit) and constraints like "coefficients sum to 1" in constrained portfolios.

### Part (b): Inequality constraints $A\beta \leq b$ (or $\geq b$) — NO general closed form

With inequalities the problem becomes **quadratic programming (QP)**: minimize a convex quadratic subject to linear inequalities. There is **no general closed-form solution** — it must be solved numerically (active-set, interior-point, or projected-gradient methods). The reason is structural:

**Summary table:**

| Constraint | Solvable? | Form |
|---|---|---|
| None (OLS) | Yes | $\hat\beta = (X^\top X)^{-1}X^\top y$ |
| $A\beta = b$ (equality) | Yes — **closed form** | Projection formula (boxed above), via Lagrange multipliers |
| $A\beta \leq b$ / $\geq b$ (inequality) | Yes, but **no closed form** | Quadratic program; closed-form *only conditional on the active set*, found numerically via KKT |

---

## Problem 5: If you duplicate the dataset m times, how do beta, t-stats, p-values, and R-squared change?

**Problem:** You run OLS of $y$ on $X$ with $n$ observations. Now you **stack the exact same dataset $m$ times** (so the new sample is $mn$ rows, each original row repeated $m$ times) and re-run OLS. How do $\hat\beta$, the t-statistics, the p-values, and $R^2$ change?

**Answer summary:** $\hat\beta$ and $R^2$ are **unchanged**; standard errors **shrink by $\sqrt{m}$**, so t-stats **inflate by $\sqrt{m}$** and p-values **collapse toward zero**. This is the classic "fake sample size" trap — duplicating data fabricates statistical significance without adding any real information.

**Setup:** Stacking $m$ copies gives $\tilde X = \mathbf{1}_m \otimes X$ and $\tilde y = \mathbf{1}_m \otimes y$. The cross-products simply scale by $m$:

$$\tilde X^\top \tilde X = m\,X^\top X, \qquad \tilde X^\top \tilde y = m\,X^\top y$$

### $\hat\beta$ — unchanged

$$\hat\beta_{\text{new}} = (\tilde X^\top \tilde X)^{-1}\tilde X^\top \tilde y = (m X^\top X)^{-1}(m X^\top y) = (X^\top X)^{-1}X^\top y = \hat\beta_{\text{old}}$$

The factor $m$ cancels. Intuitively, the *minimizer* of the squared error doesn't move — scaling every squared residual by $m$ scales the whole objective but not its argmin. Each residual $\hat\varepsilon_i$ is identical to before (just repeated $m$ times).

### $R^2$ — unchanged

$\text{SS}_{\text{res}}$ and $\text{SS}_{\text{tot}}$ both scale by exactly $m$ (each squared term is repeated $m$ times), so their ratio is invariant:

$$R^2_{\text{new}} = 1 - \frac{m\,\text{SS}_{\text{res}}}{m\,\text{SS}_{\text{tot}}} = R^2_{\text{old}}$$

(Adjusted $R^2$ shifts microscopically because its $n$ changes, but ordinary $R^2$ is identical.)

### Standard errors — shrink by $\sqrt{m}$

This is where the damage is. Two things change in the SE formula $\text{SE}(\hat\beta_j) = \sqrt{\hat\sigma^2 \,[(X^\top X)^{-1}]_{jj}}$:

1. **The $(X^\top X)^{-1}$ factor shrinks:** $(\tilde X^\top \tilde X)^{-1} = \frac{1}{m}(X^\top X)^{-1}$.
2. **The variance estimate $\hat\sigma^2$ barely changes** — and this is the crux. The naive estimator uses the *inflated* degrees of freedom:

$$\hat\sigma^2_{\text{new}} = \frac{\text{SS}_{\text{res,new}}}{mn - p} = \frac{m\,\text{SS}_{\text{res,old}}}{mn - p} \;\approx\; \frac{m\,\text{SS}_{\text{res,old}}}{mn} = \frac{\text{SS}_{\text{res,old}}}{n} \approx \hat\sigma^2_{\text{old}}$$

The numerator grew by $m$ but the denominator (d.o.f.) also grew by $\approx m$, so $\hat\sigma^2$ is essentially unchanged. Combining:

$$\text{Var}(\hat\beta_{\text{new}}) = \hat\sigma^2_{\text{new}}\,(\tilde X^\top\tilde X)^{-1} \approx \hat\sigma^2_{\text{old}}\cdot \tfrac{1}{m}(X^\top X)^{-1} = \frac{1}{m}\text{Var}(\hat\beta_{\text{old}})$$

$$\boxed{\text{SE}_{\text{new}} \approx \frac{\text{SE}_{\text{old}}}{\sqrt{m}}}$$

### t-stats — inflate by $\sqrt{m}$, p-values — collapse

Since $\hat\beta$ is unchanged and the SE shrank by $\sqrt m$:

$$t_{\text{new}} = \frac{\hat\beta_j}{\text{SE}_{\text{new}}} \approx \sqrt{m}\cdot t_{\text{old}}$$

The t-stat blows up by $\sqrt m$. The **F-stat** inflates the same way. So you can manufacture "significance" for *any* nonzero coefficient just by copy-pasting rows.

**How the p-value is computed, and why it → 0.** For the two-sided test of $H_0: \beta_j = 0$, the p-value is the tail probability of the reference $t$-distribution beyond the observed statistic:

$$p =,\Pr\!\big(T_\nu > |t|\big), \qquad \nu = \text{d.o.f.} = mn - p$$

Two forces both push this tail probability toward 0 under duplication:

1. **The statistic moves out into the tail.** $|t|$ grows by $\sqrt m$, and the $t$ (and normal) tail $\Pr(T_\nu > |t|)$ decays **faster than exponentially** in $|t|$ (roughly $\sim e^{-t^2/2}$ for the normal limit). So scaling $|t|$ by $\sqrt m$ shrinks the tail probability extremely fast — e.g. a borderline $t = 2$ ($p \approx 0.046$) becomes $t = 2\sqrt{4} = 4$ at $m=4$, giving $p \approx 6\times 10^{-5}$, and $t=2\sqrt{9}=6$ at $m=9$ gives $p \approx 2\times10^{-9}$.

The dominant term is the $e^{-t^2/2}$ tail decay evaluated at $t \approx \sqrt m\, t_{\text{old}}$, so the p-value falls **roughly like $e^{-m\, t_{\text{old}}^2/2}$ — exponentially in the number of copies $m$.** That is why even a weak, barely-significant effect becomes "astronomically significant" after a handful of duplications, despite zero new information.

**Why this is a trap — the i.i.d. violation:** the math is correct but the *inference* is invalid. The SE formula assumes $mn$ **independent** observations. Duplicated rows are perfectly correlated, not independent, so the true information content is still just $n$ points. The standard errors are wrong because the independence assumption (assumption 5, no error correlation — see [[linear_regression]]) is grossly violated: $\text{Cov}(\varepsilon_i, \varepsilon_j) = \sigma^2$ for duplicate pairs, not 0.

**Connections:**
- This is the degenerate extreme of **clustered data** — each set of $m$ duplicates is a cluster with perfect within-cluster correlation. **Cluster-robust standard errors** (clustering on the original row id) would correctly recover $\text{SE}_{\text{old}}$ and undo the inflation.
- It's the mirror image of [[#Problem 3: How many points are needed to halve the standard deviation of beta of a simple linear regression?]]: there, $4\times$ *genuinely independent* points halve the SE; here, $4\times$ *duplicated* points appear to halve the SE but add zero real information.
- Practical takeaways: never let train/test or resampling pipelines duplicate rows; watch for it when joining tables (a one-to-many join silently replicates rows and fabricates significance).
