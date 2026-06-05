#tags: factor-model, IC, optimization, Lagrangian, active-management

# Combining Factors: IC Maximization

## Problem Setup

Given two factors F₁ and F₂ with:
- Std(F₁) = Std(F₂) = 1
- Cov(F₁, F₂) = 0 (orthogonal)
- IC of Fᵢ w.r.t. return r: ρᵢ = Corr(Fᵢ, r)

Form a combined factor F = w₁F₁ + w₂F₂. We want to choose w₁, w₂ to **maximize the IC of F**, subject to **Std(F) = 1**.

---

## Translating the Constraint

Since F₁, F₂ are orthogonal with unit variance:

$$\text{Var}(F) = w_1^2 \text{Var}(F_1) + w_2^2 \text{Var}(F_2) = w_1^2 + w_2^2$$

So Std(F) = 1 becomes:

$$w_1^2 + w_2^2 = 1$$

---

## Translating the Objective

$$\text{IC} = \text{Corr}(F, r) = \frac{\text{Cov}(F, r)}{\text{Std}(F) \cdot \text{Std}(r)}$$

Since Std(F) = 1 by constraint:

$$\text{IC} = \frac{\text{Cov}(F, r)}{\text{Std}(r)} = \frac{w_1 \text{Cov}(F_1, r) + w_2 \text{Cov}(F_2, r)}{\text{Std}(r)}$$

Using $\rho_i = \text{Corr}(F_i, r) = \frac{\text{Cov}(F_i, r)}{\text{Std}(r)}$:

$$\text{IC} = w_1 \rho_1 + w_2 \rho_2$$

---

## Optimization Problem

$$\max_{w_1, w_2} \quad w_1 \rho_1 + w_2 \rho_2 \qquad \text{s.t.} \quad w_1^2 + w_2^2 = 1$$

---

## Solution via Lagrangian

Form the Lagrangian:

$$\mathcal{L} = w_1 \rho_1 + w_2 \rho_2 - \lambda(w_1^2 + w_2^2 - 1)$$

First-order conditions:

$$\frac{\partial \mathcal{L}}{\partial w_1} = \rho_1 - 2\lambda w_1 = 0 \implies w_1 = \frac{\rho_1}{2\lambda}$$

$$\frac{\partial \mathcal{L}}{\partial w_2} = \rho_2 - 2\lambda w_2 = 0 \implies w_2 = \frac{\rho_2}{2\lambda}$$

Substitute into the constraint:

$$\frac{\rho_1^2 + \rho_2^2}{4\lambda^2} = 1 \implies \lambda = \pm\frac{\sqrt{\rho_1^2 + \rho_2^2}}{2}$$

Take the positive root (to maximize, not minimize):

$$\boxed{w_i^* = \frac{\rho_i}{\sqrt{\rho_1^2 + \rho_2^2}}}$$

---

## Combined IC

$$\text{IC}^* = w_1^* \rho_1 + w_2^* \rho_2 = \frac{\rho_1^2 + \rho_2^2}{\sqrt{\rho_1^2 + \rho_2^2}} = \boxed{\sqrt{\rho_1^2 + \rho_2^2}}$$

---

## Key Takeaways

1. **Optimal weights** simply normalize the IC vector — w* is the unit vector in the direction of (ρ₁, ρ₂).
2. **Combined IC combines in quadrature**: independent signals stack as root-sum-of-squares, not by simple addition.
3. **Cauchy-Schwarz interpretation**: the bound **w**·**ρ** ≤ ‖**w**‖₂‖**ρ**‖₂ = ‖**ρ**‖₂ is achieved exactly when **w** ∝ **ρ**.
4. **Fundamental Law of Active Management**: this is the multi-signal generalization — orthogonal signals with ICs ρ₁, ρ₂ yield a combined IC of √(ρ₁² + ρ₂²), which is always ≥ max(ρ₁, ρ₂).

## See Also

- [[fundamental_law_of_active_management]]
- [[factor_model]]
- [[information_coefficient]]
