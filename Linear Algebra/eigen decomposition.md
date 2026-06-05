

**eigen decomposition** is a method of breaking down a complex matrix into its fundamental, structural components.

If you think of a square matrix as a mathematical instruction that stretches, skews, and rotates space, eigen decomposition finds the specific, hidden axes in that space that _only stretch_ and _do not rotate_.

### The General Concept

For a given square matrix $A$, an eigenvector $\mathbf{v}$ and its corresponding eigenvalue $\lambda$ satisfy this elegant equation:

$$A \mathbf{v} = \lambda \mathbf{v}$$

This means that when you apply the matrix transformation $A$ to the vector $\mathbf{v}$, the vector doesn't change its direction. It simply scales by a factor of $\lambda$.

When you find all the eigenvectors and eigenvalues for the matrix, you can decompose the entire matrix into:

$$A = V \Lambda V^{-1}$$

Where $V$ is a matrix containing all the eigenvectors (each column of v is an eigen vector) , and $\Lambda$ (Lambda) is a diagonal matrix containing all the eigenvalues. You have effectively isolated the "directions" of the matrix ($V$) from the "magnitude" or scale of those directions ($\Lambda$).

### The Connection to the Covariance Matrix

The relationship between eigen decomposition and the covariance matrix, denoted as $\Sigma$, is the mathematical engine behind almost all modern quantitative risk and factor modeling.

The covariance matrix is mathematically special: it is always **symmetric** and **positive semi-definite**. Because of this, a rule called the Spectral Theorem applies, which makes its eigendecomposition incredibly clean:

$$\Sigma = V \Lambda V^T$$

The inverse $V^{-1}$ simply becomes the transpose $V^T$. More importantly, this guarantees that all the resulting eigenvectors are **perfectly orthogonal (perpendicular)** to each other.

Here is what that translates to in the real world:

- **The Eigenvectors ($V$):** These represent the independent "directions" of maximum variance in your data. Because they are orthogonal, they are statistically uncorrelated. In finance, if you run this on an asset returns covariance matrix, these eigenvectors become your purely uncorrelated, synthetic risk factors (often called Principal Components).
    
- **The Eigenvalues ($\Lambda$):** These measure the magnitude of variance that exists exactly along those eigenvector directions. A massive eigenvalue means that specific eigenvector explains a huge portion of the system's total risk (like a broad market beta factor). A tiny eigenvalue means very little risk is driven by that direction.