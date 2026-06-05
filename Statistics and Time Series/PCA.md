#tags: PCA, eigen decomposition, factor models, covariance matrix, dimensionality reduction

# Principal Component Analysis (PCA)

## Introduction

PCA is a dimensionality reduction technique that finds the directions of maximum variance in data. It rotates the original coordinate system into a new basis where axes are ordered by how much variance they explain — the first principal component captures the most, the second the next most (orthogonal to the first), and so on.

**Core idea:** Replace $p$ correlated variables with $k \ll p$ uncorrelated components that retain most of the information.
For simplicity, we can always assume that its demeaned such that var(X) = $X^TX/n$



**Steps:**
1. Compute sample covariance matrix $\hat{\Sigma}$
2. Eigen decompose: $\hat{\Sigma} = Q \Lambda Q^\top$
3. Project data onto top $k$ eigenvectors: $Z = Q_k^\top X$

Alternatively, apply **SVD** directly to the data matrix (more numerically stable):
$$X = USV^\top \implies \text{PCs} = V, \quad \lambda_i = s_i^2/(n-1)$$

---

## PCA via Eigen Decomposition

The covariance matrix $\Sigma$ decomposes as:

$$\Sigma = Q \Lambda Q^\top$$

where:
- $Q = [q_1, q_2, \ldots, q_p]$ — orthonormal eigenvectors (principal directions)
- $\Lambda = \text{diag}(\lambda_1, \lambda_2, \ldots, \lambda_p)$ — eigenvalues with $\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_p \geq 0$

**Each eigenvector** $q_i$ is a principal component direction — the axis along which the $i$-th most variance lies.

**Each eigenvalue** $\lambda_i$ is the variance of the data projected onto $q_i$.

**Variance explained by top $k$ components:**
$$\text{Explained variance} = \frac{\sum_{i=1}^k \lambda_i}{\sum_{i=1}^p \lambda_i} = \frac{\sum_{i=1}^k \lambda_i}{\text{tr}(\Sigma)}$$

**Key properties:**
$$\text{tr}(\Sigma) = \sum_i \lambda_i \quad \text{(total variance preserved)}$$
$$\det(\Sigma) = \prod_i \lambda_i$$

If any $\lambda_i = 0$, $\Sigma$ is singular — multicollinearity in the data.

---

## Factor Models & Risk Decomposition

In quant research, eigen decomposition of $\Sigma$ is the backbone of **statistical factor models**.

$$\Sigma \approx Q_k \Lambda_k Q_k^\top + D$$

where:
- $Q_k \Lambda_k Q_k^\top$ — systematic risk explained by top $k$ factors
- $D$ — diagonal matrix of idiosyncratic (stock-specific) variances

**Each eigenvector = a latent risk factor.** The top eigenvectors typically correspond to interpretable economic factors:

**Eigenvalues tell you how much risk each factor carries:**

$$\text{Risk contribution of factor } i = \frac{\lambda_i}{\sum_j \lambda_j}$$

**Portfolio variance decomposition:**

$$\text{Var}(w^\top r) = w^\top \Sigma w = \sum_{i=1}^k \lambda_i (q_i^\top w)^2 + w^\top D w$$

The term $(q_i^\top w)^2$ is the portfolio's **loading** on factor $i$ — how much exposure the portfolio has to that risk factor.

**Practical uses at a firm like Acadian:**
- Identify dominant risk factors driving cross-sectional return variation
- Construct market-neutral or factor-neutral portfolios by zeroing out loadings on top PCs
- Regularize noisy sample covariance matrices by truncating small eigenvalues (eigenvalue clipping) or applying Ledoit-Wolf shrinkage

---

## Related Notes
- [[variance_covariance_correlation]]
- [[linear_regression]]
