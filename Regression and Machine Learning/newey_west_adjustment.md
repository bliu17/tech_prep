#tags: regression, newey-west, HAC, autocorrelation, heteroskedasticity, standard-errors, time-series

# Newey-West (HAC) Adjustment

Matrix-form construction of the heteroskedasticity- and autocorrelation-consistent (HAC) standard errors. This is the detail behind assumption 5.4 in [[linear_regression]].

### 1. The Core Target

Every standard error calculation relies on the "sandwich" estimator for the variance of the coefficients:

$$Var(\hat{\beta}) = (X^T X)^{-1} S (X^T X)^{-1}$$

The only thing that changes between OLS, White (Robust), and Newey-West is how you construct the $K \times K$ "meat" matrix, $S = X^T \Omega X$.

### 2. The Weight Matrix ($W$)

Instead of summations and lags, Newey-West defines a $T \times T$ symmetric, banded Toeplitz matrix, $W$. This matrix maps the temporal "memory" of the data based on your chosen lag limit, $L$.

- The main diagonal is exactly $1$.
    
- The off-diagonals decay linearly: $W_{i,j} = 1 - \frac{|i-j|}{L+1}$ for all $|i-j| \le L$.
    
- Every entry where $|i-j| > L$ is exactly $0$.
    

### 3. The Residual Matrix ($E$)

Take your $T \times 1$ vector of standard OLS residuals ($\hat{e}$) and project them onto the diagonal of a $T \times T$ matrix, $E$:

$$E = diag(\hat{e}_1, \hat{e}_2, \dots, \hat{e}_T)$$

### 4. The Newey-West Implied Covariance ($\Omega_{NW}$)

Using pure matrix multiplication, the Newey-West implied $T \times T$ residual covariance matrix is constructed by "sandwiching" the temporal weights between the residuals:

$$\hat{\Omega}_{NW} = E W E$$

_Notice how this works geometrically:_ Multiplying $E \times W \times E$ scales every temporal weight $W_{i,j}$ by the outer product of the residuals at those specific times ($\hat{e}_i \hat{e}_j$).

### 5. The Final Correction

Now, plug that implied covariance matrix back into the meat of the original equation.

The fully corrected Newey-West variance-covariance matrix of your factor estimates is:

$$Var_{NW}(\hat{\beta}) = (X^T X)^{-1} [X^T (E W E) X] (X^T X)^{-1}$$

**Why this is mathematically brilliant:**

If you set $L=0$ (assuming no autocorrelation), the weight matrix $W$ collapses into an Identity matrix ($I$). $E I E$ simply becomes $E^2$, which gives you the standard White (HC0) heteroskedasticity-robust standard errors. If you assume constant variance on top of that, it collapses all the way back down to the classic OLS formula $\sigma^2(X^TX)^{-1}$.


### 6. Example of W when lag = 2 and 3

lag = 2
$$W = \begin{bmatrix} 1.00 & 0.67 & 0.33 & 0 & 0 & 0 & \dots \\ 0.67 & 1.00 & 0.67 & 0.33 & 0 & 0 & \dots \\ 0.33 & 0.67 & 1.00 & 0.67 & 0.33 & 0 & \dots \\ 0 & 0.33 & 0.67 & 1.00 & 0.67 & 0.33 & \dots \\ 0 & 0 & 0.33 & 0.67 & 1.00 & 0.67 & \dots \\ 0 & 0 & 0 & 0.33 & 0.67 & 1.00 & \dots \\ \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \ddots \end{bmatrix}$$

lag = 3
$$W = \begin{bmatrix} 1.00 & 0.75 & 0.50 & 0.25 & 0 & 0 & \dots \\ 0.75 & 1.00 & 0.75 & 0.50 & 0.25 & 0 & \dots \\ 0.50 & 0.75 & 1.00 & 0.75 & 0.50 & 0.25 & \dots \\ 0.25 & 0.50 & 0.75 & 1.00 & 0.75 & 0.50 & \dots \\ 0 & 0.25 & 0.50 & 0.75 & 1.00 & 0.75 & \dots \\ 0 & 0 & 0.25 & 0.50 & 0.75 & 1.00 & \dots \\ \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \ddots \end{bmatrix}$$

## Related Notes

- [[linear_regression]] — assumption 5 (no autocorrelation), where this correction is applied
- [[GLS]] — the general framework; White, WLS, and Fama-MacBeth as special cases