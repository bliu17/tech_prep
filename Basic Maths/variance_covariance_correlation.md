#tags: variance, covariance, correlation, statistics, probability

# Variance, Covariance & Correlation

## Variance

$$\text{Var}(X) = E[X^2] - (E[X])^2$$

$$\text{Var}(X) = E[(X - \mu)^2]$$

**Sample variance (unbiased):**
$$s^2 = \frac{1}{n-1}\sum_{i=1}^n (x_i - \bar{x})^2$$

**Linear transformation:**
$$\text{Var}(aX + b) = a^2 \text{Var}(X)$$

**Sum of two variables:**
$$\text{Var}(X + Y) = \text{Var}(X) + \text{Var}(Y) + 2\text{Cov}(X, Y)$$

**Sum of $n$ independent variables:**
$$\text{Var}\!\left(\sum_{i=1}^n X_i\right) = \sum_{i=1}^n \text{Var}(X_i) \quad \text{(if independent)}$$

---

## Covariance

$$\text{Cov}(X, Y) = E[XY] - E[X]E[Y]$$

$$\text{Cov}(X, Y) = E[(X - \mu_X)(Y - \mu_Y)]$$

**Sample covariance:**
$$\hat{\sigma}_{XY} = \frac{1}{n-1}\sum_{i=1}^n (x_i - \bar{x})(y_i - \bar{y})$$

**Key properties:**

$$\text{Cov}(X, X) = \text{Var}(X)$$

$$\text{Cov}(aX, bY) = ab \cdot \text{Cov}(X, Y)$$

$$\text{Cov}(X + Y, Z) = \text{Cov}(X, Z) + \text{Cov}(Y, Z)$$

---

## Correlation

$$\rho(X, Y) = \frac{\text{Cov}(X, Y)}{\sigma_X \sigma_Y}, \quad \rho \in [-1, 1]$$

**Sample correlation:**
$$\hat{\rho} = \frac{\hat{\sigma}_{XY}}{\hat{\sigma}_X \hat{\sigma}_Y}$$

**Bound (from Cauchy-Schwarz):**
$$|\text{Cov}(X,Y)|^2 \leq \text{Var}(X)\cdot\text{Var}(Y)$$

---

## Portfolio / Matrix Form

**Portfolio variance:**
$$\text{Var}(w^\top r) = w^\top \Sigma w$$

**Covariance matrix:**
$$\Sigma_{ij} = \text{Cov}(r_i, r_j)$$

**$\Sigma$ must be positive semi-definite:**
$$w^\top \Sigma w \geq 0 \quad \forall\, w$$

---

## Key Facts
Cov$(X,Y) = 0$ $\Rightarrow$ Independent ❌
Independent $\Rightarrow$ Cov$(X,Y) = 0$✅
**Counterexample for zero cov ≠ independence:** Let $X \sim \text{Uniform}(-1,1)$, $Y = X^2$. Then $\text{Cov}(X,Y) = 0$ but $Y$ is fully determined by $X$.

**Proof:**
$$\text{Cov}(X, Y) = E[XY] - E[X]E[Y]$$
- $E[X] = 0$
- $E[XY] = E[X \cdot X^2] = E[X^3] = \int_{-1}^{1} x^3 \cdot \frac{1}{2}\,dx = 0$ 
$$\therefore \text{Cov}(X, Y) = 0 - 0 = 0$$

---

## Matrix Forms

*Throughout, assume the data is demeaned: $X$ is an $n \times p$ data matrix ($n$ observations, $p$ variables) whose columns have zero mean, and $x$ is a $p \times 1$ random vector with $E[x] = 0$.*

**Core definitions:**
$$\Sigma = \frac{1}{n}X^\top X$$

Entry $(j,k)$ is $\text{Cov}(\text{col}_j, \text{col}_k)$; the diagonal holds the individual variances. Use $\tfrac{1}{n-1}$ for the unbiased estimator.

**Linear transformation (the workhorse):**
$$\text{Cov}(Ax) = A\,\Sigma\,A^\top, \qquad \text{Cov}(Ax,\, By) = A\,\text{Cov}(x,y)\,B^\top$$

Scalar special case $\text{Var}(w^\top x) = w^\top \Sigma w$ is just $A = w^\top$.

**Portfolio / factor algebra:**
$$\sigma_p^2 = w^\top \Sigma w, \qquad \text{Cov}(w_1^\top x,\, w_2^\top x) = w_1^\top \Sigma w_2$$

Factor model $x = Bf + \varepsilon$ ($B$ = loadings, $f$ = factors, $\varepsilon \perp f$):
$$\Sigma = B\,\Sigma_f\,B^\top + D \quad (D \text{ diagonal idiosyncratic variances})$$

**Decompositions:**
$$\Sigma = Q\Lambda Q^\top \quad (\text{eigendecomp, } Q \text{ orthonormal}), \qquad \Sigma = LL^\top \quad (\text{Cholesky})$$

See [[eigen decomposition]] for the full treatment.

Cholesky generates correlated normals: if $z \sim N(0, I)$ then $Lz \sim N(0, \Sigma)$. The precision matrix $\Sigma^{-1}$ appears in the Gaussian density and GLS; its off-diagonals encode partial (conditional) correlations.

**Correlation in matrix form:**
$$P = D_\sigma^{-1/2}\,\Sigma\,D_\sigma^{-1/2}, \qquad D_\sigma = \text{diag}(\Sigma)$$

**Dual view ($X^\top X$ vs $X X^\top$):**
- $X^\top X$ is $p \times p$ — covariance between **variables** (the matrix for $w^\top \Sigma w$).
- $X X^\top$ is $n \times n$ — Gram matrix, covariance between **observations** (sample-space / kernel view).
- Both are symmetric PSD and share the **same nonzero eigenvalues**. Via SVD $X = U\Sigma_{\text{svd}} V^\top$: columns of $V$ are eigenvectors of $X^\top X$, columns of $U$ are eigenvectors of $X X^\top$. Compute PCA from whichever is smaller.

**Properties to state cold:**
- $\Sigma$ symmetric PSD: $w^\top \Sigma w \geq 0 \;\; \forall w$; strictly PD iff full-rank (no perfect collinearity).
- $\text{tr}(\Sigma) = \sum_i \lambda_i$ = total variance (the quantity PCA partitions).
- $\text{rank}(\Sigma) = \text{rank}(X)$.

---

