#tags: regression, ridge, lasso, regularization, shrinkage, linear-algebra

# Ridge and Lasso Regression

Both methods add a penalty to the OLS objective to shrink coefficients — addressing overfitting and multicollinearity. Both always introduce bias in the $\hat{\beta}$ estimates, ie. if $\beta>0$ then the bias drives the new $\beta$ estimate to be smaller than true beta, and vice versa.

---

## Ridge Regression (L2)

**Objective:**

$$\min_\beta \|y - X\beta\|^2 + \lambda \|\beta\|^2 = \min_\beta (y - X\beta)^\top(y - X\beta) + \lambda \beta^\top \beta$$

### Closed-Form Solution

Taking the derivative and setting to zero:

$$\frac{\partial}{\partial \beta}\left[(y - X\beta)^\top(y - X\beta) + \lambda \beta^\top \beta\right] = -2X^\top(y - X\beta) + 2\lambda\beta = 0$$

$$X^\top X \hat{\beta} + \lambda \hat{\beta} = X^\top y$$

$$(X^\top X + \lambda I)\hat{\beta} = X^\top y$$

$$\boxed{\hat{\beta}_{\text{ridge}} = (X^\top X + \lambda I)^{-1} X^\top y}$$

Compare to OLS: $\hat{\beta}_{\text{OLS}} = (X^\top X)^{-1} X^\top y$. Ridge adds $\lambda I$ to the diagonal, making the matrix always invertible regardless of multicollinearity.

### Properties

- **Biased** but lower variance than OLS — bias-variance tradeoff controlled by $\lambda$
- Shrinks all coefficients toward zero but never exactly to zero

---

## Lasso Regression (L1)

**Objective:**

$$\min_\beta \|y - X\beta\|^2 + \lambda \|\beta\|_1 = \min_\beta \sum_i(y_i - x_i^\top\beta)^2 + \lambda \sum_j |\beta_j|$$

### No Closed-Form Solution

The L1 penalty is not differentiable at $\beta_j = 0$, so there is no closed-form solution. Solved numerically via coordinate descent or LARS.

### Properties

- **Sparse solutions** — can shrink coefficients exactly to zero, performing automatic variable selection

---

## Ridge vs Lasso

|             | Ridge (L2)                            | Lasso (L1)                            |
| ----------- | ------------------------------------- | ------------------------------------- |
| Penalty     | $\lambda\|\beta\|^2$                  | $\lambda\|\beta\|_1$                  |
| Closed form | Yes                                   | No                                    |
| Sparsity    | No — shrinks but never zeros          | Yes — exact zeros                     |
| Best for    | Many small effects, multicollinearity | Few large effects, variable selection |
| Geometry    | Circular constraint region            | Diamond constraint region             |
|             |                                       |                                       |

---

## Related Notes

- [[linear_regression]]
