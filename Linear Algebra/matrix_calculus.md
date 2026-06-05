#tags: linear-algebra, calculus, gradients, matrix-calculus, optimization

---

## The Differential Method

The cleanest way to compute matrix gradients — no index juggling required.

**Goal:** derive $\partial f / \partial \beta$ where $f$ is a scalar function of vector $\beta$.

**Method:**
1. Compute the scalar differential $df$ by perturbing $\beta \to d\beta$
2. Massage $df$ into the form $g^T d\beta$
3. Read gradient directly: $\partial f / \partial \beta = g$

**Why it works:** $df$ is the first-order change in $f$. Since $f$ is scalar, $df$ must be a scalar too — so it must look like an inner product $\langle g, d\beta \rangle = g^T d\beta$. That $g$ is the gradient by definition.

---

## Another Trick: "A scalar equals its own transpose"

If $s$ is a scalar, then $s = s^T$. This lets you flip terms to collect $d\beta$ on the right:

$$(d\beta)^T A \beta = \left[(d\beta)^T A \beta\right]^T = \beta^T A^T d\beta$$

This is the workhorse step in almost every derivation.

---

## Differential Rules

| Expression | Differential |
|---|---|
| $A\beta$ | $A\, d\beta$ |
| $\beta^T A$ | $(d\beta)^T A$ |
| $\beta^T A \beta$ | $\beta^T(A + A^T)\, d\beta \quad \to 2\beta^T A\, d\beta$ if $A$ symmetric |
| $\text{tr}(A)$ | $\text{tr}(dA)$ |
| $\log \det A$ | $\text{tr}(A^{-1} dA)$ |

---

## Worked Examples

### Example 1: Quadratic form $f = \beta^T \Sigma \beta$

$$df = (d\beta)^T \Sigma \beta + \beta^T \Sigma\, d\beta$$

Apply the scalar-transpose trick to the first term:

$$(d\beta)^T \Sigma \beta = \beta^T \Sigma^T d\beta$$

So:

$$df = \beta^T \Sigma^T d\beta + \beta^T \Sigma\, d\beta = \beta^T(\Sigma + \Sigma^T)\, d\beta$$

If $\Sigma$ is symmetric ($\Sigma = \Sigma^T$):

$$df = 2\beta^T \Sigma\, d\beta \implies \boxed{\nabla_\beta f = 2\Sigma\beta}$$

**Interpretation:** portfolio variance $\beta^T \Sigma \beta$ has gradient $2\Sigma\beta$ — the marginal risk contribution of each asset.

---

### Example 2: OLS loss $f = \|y - X\beta\|^2 = (y - X\beta)^T(y - X\beta)$

Let $r = y - X\beta$, so $dr = -X\, d\beta$.

$$df = 2r^T\, dr = 2(y - X\beta)^T(-X\, d\beta) = -2(y - X\beta)^T X\, d\beta$$

So:

$$\nabla_\beta f = -2X^T(y - X\beta)$$

Setting to zero: $X^T X \hat\beta = X^T y \implies \boxed{\hat\beta = (X^T X)^{-1} X^T y}$

---

### Example 3: Linear function $f = c^T \beta$

$$df = c^T d\beta \implies \boxed{\nabla_\beta f = c}$$

---

## Matrix-to-Polynomial Conversion (Quick Reference)

See [[matrix_polynomial_forms]] for expanding $X^TX$, $X^Ty$, and quadratic forms into explicit sums.

| Matrix form | Scalar/polynomial meaning |
|---|---|
| $x^T y$ | $\sum_i x_i y_i$ |
| $x^T A x$ | $\sum_{i,j} A_{ij} x_i x_j$ |
| $(X^T X)_{ij}$ | $\sum_k X_{ki} X_{kj}$ (inner products of columns) |
| $(X^T y)_i$ | $\sum_k X_{ki} y_k$ |
| $\text{tr}(AB)$ | $\sum_{i,j} A_{ij} B_{ji}$ |

---

## Summary of derive $\partial f / \partial \beta$ which is just the gradient


1. try to find $df$ 
2. Use $s = s^T$ to flip $(d\beta)^T$ terms to $d\beta$ term 
3. derive $df = g^T d\beta$ → gradient is $g$
