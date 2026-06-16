#tags: regression, OLS, assumptions, diagnostics, beta, R-squared

# Linear Regression

## Contents

- [[#OLS Assumptions]]
- [[#OLS Estimator Formula]]
- [[#R², T-statistics, F-statistics]]
- [[#Combining OLS Estimates from Separate Datasets]]

## OLS Assumptions

The Gauss-Markov theorem guarantees OLS is BLUE (Best Linear Unbiased Estimator) when assumptions 1–5 hold. Normality (6) is needed for exact finite-sample inference.

### 1. Linearity

#### 1.1 Assumption
The true relationship is linear in parameters: $Y = X\beta + \varepsilon$.

#### 1.2 How to Check
- Scatter plots of $Y$ vs each $X_j$
- Residual vs fitted plot — should show no pattern

#### 1.3 Consequence of Violation
The model is fundamentally mis-specified

#### 1.4 What to Do if Violated
- Add polynomial terms or interaction terms
- Apply transformations (log, sqrt) to $Y$ or $X$
- Use a non-linear model

### 2. No Perfect Multicollinearity

#### 2.1 Assumption
Regressors are not perfectly linearly dependent.

#### 2.2 How to Check
- Variance Inflation Factor (VIF): VIF > 10 is a common warning threshold
- Correlation matrix of predictors

#### 2.3 Consequence of Violation
Coefficients remain unbiased but standard errors inflate, making inference unreliable. In the perfect collinearity case, $(X^\top X)$ is not invertible and OLS has no unique solution.

**Why multicollinearity inflates standard errors:**

Assume $X_1, X_2$ are demeaned. The variance of $\hat{\beta}_1$ is:

$$\text{Var}(\hat{\beta}_1) = \frac{\sigma^2}{\sum x_{i1}^2 \cdot (1 - \rho_{12}^2)}$$

where $\rho_{12} = \text{Corr}(X_1, X_2)$. As $\rho_{12} \to \pm 1$, the term $(1 - \rho_{12}^2) \to 0$, so $\text{Var}(\hat{\beta}_1) \to \infty$. 

This is exactly what VIF captures: $\text{VIF}_1 = \frac{1}{1 - R_1^2}$, where $R_1^2$ is the R² from regressing $X_1$ on all other regressors, so $\text{Var}(\hat{\beta}_1) = \frac{\sigma^2}{\sum x_{i1}^2} \cdot \text{VIF}_1$.

#### 2.4 What to Do if Violated
- Drop one of the collinear variables / Orthogonalize one to other variables, see [[Frisch-Waugh-Lovell]]
- Ridge / Lasso regression, see [[ridge_lasso]]
- PCA, see [[PCA]]

### 3. Exogeneity (Zero Conditional Mean)

#### 3.1 Assumption
**Strict form:** $E[\varepsilon \mid X] = 0$ — required for OLS to be unbiased.

**Weak form:** $\text{Cov}(\varepsilon, X) = 0$ — sufficient for OLS to be consistent, i.e. $\hat{\beta} \to \beta$ as $n \to \infty$. The strict form implies the weak form, but not vice versa.

#### 3.2 How to Check
- Hard to test directly
- Check for omitted variables, if any

#### 3.3 Consequence of Violation
Coefficients are biased and inconsistent — the most serious violation as it cannot be fixed by larger samples. To solve questions towards this, we can assume the true model being $y = X\beta +Z\gamma +u$ where u is the true residual. since $\epsilon = Z\gamma +u$, we can derive that $E[\beta] = \beta + (X^TX)^{-1}X^T Z\gamma$ and if we think X and Z being demeaned, $Cov(X,Z) = X^T Z$ >0 $\rightarrow E[\epsilon|X] >0$, beta will bias up.

#### 3.4 What to Do if Violated
- Find the missing variables if available. This one is hard to test and fix .....

### 4. Homoskedasticity

#### 4.1 Assumption
$\text{Var}(\varepsilon_i \mid X) = \sigma^2$ — constant error variance across observations.

#### 4.2 How to Check
- Residual vs fitted plot — look for a funnel shape
- Breusch-Pagan test
- White test

#### 4.3 Consequence of Violation
Coefficients remain unbiased but OLS is no longer efficient. Standard errors are incorrect (typically understated), so t-stats and p-values are unreliable. True Standard errors of beta are larger than estimates when the extreme data points (in terms of X) tends to have larger error, which is often the case. True standard errors of beta are smaller while the points near center mu(x) have larger error, which is rare.

#### 4.4 What to Do if Violated
- Use **heteroskedasticity-robust (White) standard errors**
- WLS (Weighted Least Squares) , weight proportion to residual variance, ie. $w_i  = 1 / \sigma^2_i$

**White robust standard error

$$\text{Var}(\hat{\beta} \mid X) = (X^\top X)^{-1} \, X^\top \Omega X \, (X^\top X)^{-1}.$$

 $\Omega$ is unknown, and it has $n$ distinct diagonal entries $\sigma_i^2$ (one per observation).We estimate each $\sigma_i^2$ with the squared OLS residual $\hat{\varepsilon}_i^2$:

$$\widehat{\text{Var}}(\hat{\beta})_{HC0} = (X^\top X)^{-1} \, (X^\top \hat{\Omega} X) \, (X^\top X)^{-1}, \qquad X^\top \hat{\Omega} X = \sum_{i=1}^n \hat{\varepsilon}_i^2 \, x_i x_i^\top.$$

A single $\hat{\varepsilon}_i^2$ is a terrible estimate of $\sigma_i^2$ (it's based on one point), but the *sum* $\sum_i \hat{\varepsilon}_i^2 x_i x_i^\top$ converges to $\sum_i \sigma_i^2 x_i x_i^\top$ by a law of large numbers.
**Finite-sample variants.** Because OLS residuals are systematically too small (they're fit to the data), HC0 understates variance in small samples. Common corrections rescale the residuals:

| Estimator | Meat weight on $\hat{\varepsilon}_i^2$   | Note                                                              |
| --------- | ---------------------------------------- | ----------------------------------------------------------------- |
| HC0       | $\hat{\varepsilon}_i^2$                  | White's original; biased down in small $n$                        |
| HC1       | $\frac{n}{n-p}\,\hat{\varepsilon}_i^2$   | Simple d.o.f. scaling; Stata's default                            |
| HC2       | $\hat{\varepsilon}_i^2 / (1 - h_{ii})$   | $h_{ii}$ = leverage; exactly unbiased under homoskedasticity      |
| HC3       | $\hat{\varepsilon}_i^2 / (1 - h_{ii})^2$ | Approximates jackknife; preferred when high-leverage points exist |

where $h_{ii} = x_i^\top (X^\top X)^{-1} x_i$ is the $i$-th hat-matrix diagonal. Use **HC3** as the default for small/moderate samples; HC0/HC1 are fine when $n$ is large.

**Caveats.** Robust SEs fix *inference*, not the estimator — $\hat{\beta}$ is still unbiased but inefficient, so if you know the variance structure, WLS is better. 

### 5. No Autocorrelation

#### 5.1 Assumption
$\text{Cov}(\varepsilon_i, \varepsilon_j) = 0$ for $i \neq j$ — errors are uncorrelated across observations. Especially relevant in time series.

#### 5.2 How to Check
- Durbin-Watson statistic (DW ≈ 2 → no autocorrelation; DW < 2 → positive autocorrelation)
- Breusch-Godfrey test (more general, handles higher-order autocorrelation)
- ACF/PACF plot of residuals

#### 5.3 Consequence of Violation
Coefficients remain unbiased but OLS is inefficient. Standard errors are typically underestimated as typically we find that  $\text{Cov}(\varepsilon_i, \varepsilon_j) >0$ , leading to inflated t-stats.

#### 5.4 What to Do if Violated
- Use **Newey-West (HAC) standard errors**, use sample residuals from the OLS to calculate an estimate for the sample autocovariances for a specific number of lags. For the full matrix-form derivation of the HAC correction (the banded Toeplitz weight matrix $W$, $\hat{\Omega}_{NW}=EWE$, and how it collapses to White/OLS), see [[newey_west_adjustment]].
- GLS
- if it's panel data pool regression, we can add lagged dependent variable
- If it's panel data, consider using Fama-MacBeth regression.
- If it's panel data, consider clustered standard error.

For the **Newey-West (HAC)** estimator and the general **GLS** framework — and how White, WLS, and Fama-MacBeth are all special cases or approximations of GLS — see [[GLS and OLS with Robust Estimation]].

> **Note:** Assumptions 4 and 5 together are equivalent to the single matrix condition $\text{Var}(\varepsilon \mid X) = \sigma^2 I$ — homoskedasticity gives equal diagonal entries $\sigma^2$, and no autocorrelation gives zero off-diagonal entries.

### 6. Normality of Errors (not strictly needed)

#### 6.1 Assumption
$\varepsilon \mid X \sim N(0, \sigma^2)$. Not required for unbiasedness or consistency, but needed for exact t- and F-statistics in finite samples.

#### 6.2 How to Check
- Q-Q plot of residuals

#### 6.3 Consequence of Violation
t- and F-statistics no longer follow their exact distributions, so p-values are approximate. In large samples this is minor due to the [[central_limit_theorem|CLT]].

**Why the CLT rescues inference even when errors aren't normal.** The OLS estimator is a **linear combination of the errors**:

$$\hat{\beta} = \beta + (X^\top X)^{-1} X^\top \varepsilon = \beta + \sum_{i=1}^n w_i \, \varepsilon_i$$

i.e. a *weighted sum* of the $n$ error terms (with weights $w_i$ from the rows of $(X^\top X)^{-1}X^\top$). The CLT applies to exactly this kind of sum: as $n \to \infty$, a weighted sum of i.i.d. (or independent, non-dominant) finite-variance errors is **asymptotically normal regardless of the errors' own distribution**. So:

$$\sqrt{n}\,(\hat{\beta} - \beta) \;\xrightarrow{d}\; N\!\big(0,\; \sigma^2 \, Q^{-1}\big), \qquad Q = \lim \tfrac{1}{n}X^\top X$$

Normality of $\hat\beta$ is what the t- and F-tests actually require — not normality of the raw errors. Two consequences:

- The **t-statistic** $\hat\beta_j / \text{SE}(\hat\beta_j)$ is asymptotically standard normal, and since $t_{n-p-1} \to N(0,1)$ as $n$ grows, the usual $t$ critical values stay approximately valid.
- The **F-statistic** $\times\, p$ is asymptotically $\chi^2_p$ (a Wald statistic), so the F-test remains approximately valid too.

**Caveat (the Lindeberg condition).** This needs no single observation to dominate the weighted sum

---

### Summary Table

| # | Assumption | Key Test | Fix |
|---|-----------|----------|-----|
| 1 | Linearity | Residual vs fitted plot | Transform variables, add nonlinear terms |
| 2 | No multicollinearity | VIF | Ridge, drop variables, PCA |
| 3 | Exogeneity | Hausman test, domain knowledge | IV/2SLS, fixed effects |
| 4 | Homoskedasticity | Breusch-Pagan, White test | Robust SEs, WLS |
| 5 | No autocorrelation | Durbin-Watson, Breusch-Godfrey | Newey-West SEs, GLS |
| 6 | Normality of errors | Q-Q plot, Jarque-Bera | Bootstrap, large-sample CLT |

---

## OLS Estimator Formula

### Simple Linear Regression

For $Y = \beta_0 + \beta_1 X + \varepsilon$, OLS minimizes the sum of squared residuals:

$$\hat{\beta}_1 = \frac{\text{Cov}(X, Y)}{\text{Var}(X)}$$

$$\hat{\beta}_0 = \bar{y} - \hat{\beta}_1 \bar{x}$$

### Multivariable Linear Regression

For $Y = X\beta + \varepsilon$ where $X$ is $n \times p$, the OLS solution is:

$$\hat{\beta} = (X^\top X)^{-1} X^\top y$$

This requires $X^\top X$ to be invertible (no perfect multicollinearity).

## R², T-statistics, F-statistics

### Simple Linear Regression

$$R^2 = \text{Corr}(X, Y)^2 = \text{Corr}(Y, \hat{Y})^2 = 1 - \frac{\text{SS}_{\text{res}}}{\text{SS}_{\text{tot}}}$$

> **Why they're equal:** Since $\hat{Y} = \hat{\beta}_0 + \hat{\beta}_1 X$ is a linear transformation of $X$, correlation is invariant (except for sign change, if  $\beta_1<0$) to linear transformations, so it gives $\text{Corr}(Y, \hat{Y})^2 = \text{Corr}(Y, X)^2$.

The t-statistic for testing $H_0: \beta_1 = 0$:

$$t = \frac{\hat{\beta}_1}{\text{SE}(\hat{\beta}_1)}, \qquad \text{SE}(\hat{\beta}_1) = \sqrt{\frac{\hat{\sigma}^2}{\sum(x_i - \bar{x})^2}}$$

where $\hat{\sigma}^2 = \frac{\text{SS}_{\text{res}}}{n - 2}$ is the estimated error variance. Under $H_0$, $t \sim t_{n-2}$.

### Multivariable Linear Regression

$$R^2 = \text{Corr}(Y, \hat{Y})^2 = 1 - \frac{\text{SS}_{\text{res}}}{\text{SS}_{\text{tot}}}$$

$$\text{SS}_{\text{res}} = (y - X\hat{\beta})^\top (y - X\hat{\beta}) = y^\top y - \hat{\beta}^\top X^\top y$$

$$\text{SS}_{\text{tot}} =  \sum_i (y_i - \bar{y})^2$$

**Adjusted R²** penalizes for adding extra predictors (use when comparing models with different numbers of regressors):

$$\bar{R}^2 = 1 - \frac{\text{SS}_{\text{res}} / (n - p - 1)}{\text{SS}_{\text{tot}} / (n - 1)}$$

The t-statistic for testing $H_0: \beta_j = 0$:

$$t = \frac{\hat{\beta}_j}{\text{SE}(\hat{\beta}_j)}, \qquad \text{SE}(\hat{\beta}_j) = \sqrt{\hat{\sigma}^2 \cdot [(X^\top X)^{-1}]_{jj}}$$

where $\hat{\sigma}^2 = \frac{\text{SS}_{\text{res}}}{n - p - 1}$ and $[(X^\top X)^{-1}]_{jj}$ is the $j$-th diagonal element of $(X^\top X)^{-1}$. Under $H_0$, $t \sim t_{n-p-1}$.

Or
$$Var(\hat{\beta}|X) = \sigma^2 (X^T X)^{-1}$$

derivation:
$\hat{\beta} = (X^T X)^{-1} X^T (X\beta + \epsilon)$ -> $\hat{\beta} = \beta + (X^T X)^{-1} X^T \epsilon$
$Var(\hat{\beta}|X) = Var(\beta + (X^T X)^{-1} X^T \epsilon | X)$

$Var(\hat{\beta}|X) = Var((X^T X)^{-1} X^T \epsilon | X)$, since $Var(AZ) = A Var(Z) A^T$, 
Let $A = (X^T X)^{-1} X^T$ and let $Z = \epsilon$. Applying the rule (because $(X^T X)$ is a symmetric matrix, its inverse is also symmetric, meaning $[(X^T X)^{-1}]^T = (X^T X)^{-1}$):

$$Var(\hat{\beta}|X) = [(X^T X)^{-1} X^T] \cdot Var(\epsilon|X) \cdot [(X^T X)^{-1} X^T]^T$$
$$Var(\hat{\beta}|X) = (X^T X)^{-1} X^T \cdot Var(\epsilon|X) \cdot X (X^T X)^{-1}$$
rom our initial assumptions, we know that $Var(\epsilon|X) = \sigma^2 I$. Substitute this into the equation:

$$Var(\hat{\beta}|X) = (X^T X)^{-1} X^T (\sigma^2 I) X (X^T X)^{-1}$$

$$Var(\hat{\beta}|X) = \sigma^2 (X^T X)^{-1} X^T X (X^T X)^{-1}$$
$$Var(\hat{\beta}|X) = \sigma^2 (X^T X)^{-1}$$

and we can see that if datasets double, n goes to 2n, theoretically $(X^T X)$ will approximately double, ie. var(beta) goes to 1/2, and se(beta) goes to $1/\sqrt 2$.



**F-statistic** tests $H_0: \beta_1 = \cdots = \beta_p = 0$ jointly. Under $H_0$, $F \sim F_{p,\ n-p-1}$.

$$F = \frac{(SS_\text{tot} - SS_\text{res})/p}{SS_\text{res}/(n-p-1)} = \frac{R^2/p}{(1-R^2)/(n-p-1)} = \frac{\text{explained variance per predictor}}{\text{unexplained variance per data point with fredom}}$$

where **df** (degrees of freedom) counts the free pieces of information: the numerator uses $p$ df (one per predictor), and the denominator uses $n-p-1$ df — starting from $n$ observations and losing $p+1$ for estimating $\beta_0, \beta_1, \ldots, \beta_p$.

### Linear Regression Without Intercept

When the intercept is suppressed ($Y = \beta_1 X + \varepsilon$), the OLS estimator uses raw (uncentered) second moments:

$$\hat{\beta}_1 = \frac{\sum x_i y_i}{\sum x_i^2}$$

The standard $R^2$ can be **negative** in this case. Intuitively, forcing the line through the origin may produce a worse fit than simply predicting $\bar{y}$ for everyone — and $R^2 < 0$ is how that shows up.

---

## Combining OLS Estimates from Separate Datasets

### Setup

Two datasets with the same $p$ regressors. Stack them vertically to pool:

$$X = \begin{bmatrix} X_1 \\ X_2 \end{bmatrix}, \qquad y = \begin{bmatrix} y_1 \\ y_2 \end{bmatrix}$$

$X$ is $(n_1 + n_2) \times p$ — taller, same width. (Horizontal stacking is a join on shared observations, a different operation.)

### Pooled Estimate

Apply the standard formula $\hat{\beta} = (X^\top X)^{-1} X^\top y$. Because stacking is vertical, the cross-product sums split cleanly:

$$X^\top X = X_1^\top X_1 + X_2^\top X_2, \qquad X^\top y = X_1^\top y_1 + X_2^\top y_2$$

Note: $X^\top X_{ii}$ is the dot product of each column vectors, and the dot product is just summation of each entries square. As a result, when the matrix is stacked matrix, it will just be the summation of stacked part. Similarly for $X^\top y$.

Using each dataset's normal equations $X_j^\top y_j = (X_j^\top X_j)\hat{\beta}_j$, this becomes a **matrix-weighted average**:

$$\hat{\beta} = (A_1 + A_2)^{-1}(A_1\hat{\beta}_1 + A_2\hat{\beta}_2), \qquad A_j = X_j^\top X_j$$

You only need $A_j$ and $\hat{\beta}_j$ per dataset — never the raw rows. This is why OLS works well in distributed and streaming settings.

### Caveats

- This is an **exact algebraic identity** regardless of whether the true coefficients match across samples. It returns one pooled coefficient vector. Use a **Chow test** to check whether pooling is appropriate.
- Weighting is by $A_j = X_j^\top X_j$, **not by sample size** — a dataset with more regressor variation pulls harder. Since $\text{Var}(\hat{\beta}_j) = \sigma^2 (X_j^\top X_j)^{-1}$, this is precision-weighting, equivalent to inverse-variance weighting in meta-analysis when $\sigma^2$ is common.
- For **standard errors**, carry one extra scalar per dataset $y_j^\top y_j$. Then $\text{SSR} = (y_1^\top y_1 + y_2^\top y_2) - \hat{\beta}^\top X^\top y$, and $\hat{\sigma}^2 = \text{SSR}/(n-p)$.
- Numerically, prefer solving $(X^\top X)\hat{\beta} = X^\top y$ via Cholesky or QR over explicit inversion — same answer, better conditioning.

---

## Related Notes

- [[ridge_lasso]]
- [[fama_macbeth]]
- [[GLS and OLS with Robust Estimation]]
