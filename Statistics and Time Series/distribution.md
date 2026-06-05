#tags: distributions, normal-distribution, chi-square, statistics

# Distributions

Reference notes on common distributions used in statistical inference.

## Contents

- [Normal Distribution](#normal-distribution)
- [Chi-Square Distribution](#chi-square-distribution)
- [Student's t 分布](#students-t-分布)

---

## Normal Distribution

The **normal** (Gaussian) distribution $\mathcal{N}(\mu, \sigma^2)$ is the symmetric bell curve, fully described by its mean $\mu$ (location) and variance $\sigma^2$ (spread).

$$f(x) = \frac{1}{\sqrt{2\pi}\,\sigma}\exp\!\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)$$

Key facts:

- **Symmetric** about $\mu$; mean = median = mode.
- **Standard normal** $\mathcal{N}(0,1)$: obtained by standardizing $Z = (X-\mu)/\sigma$. CDF denoted $\Phi$.
- **68–95–99.7 rule**: ~68% of mass within $1\sigma$, ~95% within $2\sigma$ (more precisely $1.96\sigma$), ~99.7% within $3\sigma$.
- **Closed under linear combinations**: sums of independent normals are normal.
- **Why it's everywhere — the CLT**: the mean of $n$ i.i.d. finite-variance random variables tends to normal as $n\to\infty$, *regardless* of the underlying distribution. This is what justifies the z-test / t-test in practice (see [[central_limit_theorem]] and [[hypothesis_testing]]).

## Chi-Square Distribution

The **chi-square** distribution $\chi^2_k$ with $k$ degrees of freedom is the distribution of a **sum of squares of $k$ independent standard normals variable**:

$$Q = \sum_{i=1}^{k} Z_i^2 \sim \chi^2_k, \qquad Z_i \overset{\text{iid}}{\sim} \mathcal{N}(0,1)$$

Key facts:

- **Support** is $[0, \infty)$ — it's a sum of squares, so never negative.
- **Right-skewed**, especially for small $k$; becomes more symmetric (approaching normal) as $k$ grows.
- **Mean $= k$, variance $= 2k$.**
- **Degrees of freedom $k$** is the only parameter — it sets both location and shape.
- **Link to sample variance**: for i.i.d. normal data, $(n-1)s^2/\sigma^2 \sim \chi^2_{n-1}$. This is what underlies the $n-1$ degrees of freedom in the t-test and the validity of the t-statistic.
- **Uses**: variance tests, goodness-of-fit tests, tests of independence in contingency tables, and as a building block for the t-distribution (normal / $\sqrt{\chi^2/k}$) and F-distribution (ratio of two chi-squares).

## Student's t 分布

t 分布把 正态 → 卡方 → t 这条链补完整。它是用一个标准正态除以一个独立的卡方（按自由度归一化后开方）构造出来的：

$$t_\nu = \frac{Z}{\sqrt{V/\nu}}, \qquad Z \sim \mathcal{N}(0,1),\quad V \sim \chi^2_\nu,\quad Z \perp V$$

分母 $\sqrt{V/\nu}$ 是一个**随机的缩放因子**，其实就是伪装过的样本标准差（因为 $V = (n-1)s^2/\sigma^2$，所以 $\sqrt{V/\nu} = s/\sigma$）。一句话概括：

$$t = \frac{\text{标准正态}}{\text{对标准差的一个有噪声的估计（已归一化到平均值为 1）}}$$

「归一化到平均值为 1」这句话同时解释了**与正态的相似性**和**厚尾**。

### 为什么像正态（相似性）

由卡方性质，$V \sim \chi^2_\nu$ 均值 $\nu$、方差 $2\nu$，所以 $V/\nu$ 的均值 $=1$、方差 $=2/\nu$。它是一个**以 1 为中心**、离散度 $2/\nu$ 的随机量。当 $\nu$ 增大，离散度趋于 0，分母不再随机：

$$t_\nu = \frac{Z}{\sqrt{V/\nu}} \;\longrightarrow\; \frac{Z}{1} = Z \sim \mathcal{N}(0,1)$$

所以正确的心理图像是：**t 分布是一个尺度不确定的正态分布**。一旦分母稳定到平均值 1（大 $\nu$），t 就是标准正态。

### 为什么厚尾（直觉）

**直觉一——分母偶尔变小。** 分母平均是 1 但随机。$V$ 偶尔取到很小的值时，就是在用一个小于 1 的数去除 $Z$，把比值**放大**，于是极端值远比单独的 $Z$ 更常见 → 尾部更厚。关键的不对称性：除以 0.5 让数值翻倍，除以 2 只让数值减半——「小分母」事件对尾部的影响不成比例地大。除以接近零的数是爆炸性的，这就是厚尾的引擎。

**直觉二——正态的混合。** 对分母取条件后，给定某个缩放值时 $t$ 只是一个具有该尺度的正态——有时窄（分母大）、有时宽（分母小）。所以 t 是**不同宽度正态的混合**。而混合不同方差的正态，尾部总比任何单个正态都厚：宽正态成分给尾部贡献额外质量，窄正态成分把质量堆在中心，最终是**中间更尖、尾部更厚**的形状。（这也是为什么金融收益率厚尾——波动率本身每天在变，机制相同。）

### 自由度 $\nu$ 的作用

$\nu$ 控制分母晃动的程度（$\text{var}(V/\nu) = 2/\nu$）：

- **小 $\nu$** → 分母晃动剧烈 → 频繁的小分母爆炸、混合方差范围大 → **厚尾**（$\nu = 1$ 即柯西分布，厚到连均值都不存在）。
- **大 $\nu$** → 分母 ≈ 1 → 混合塌缩成单个正态 → **尾部与正态一致**。

这就是为什么 t 临界值随自由度增大而下降（df 10 时 2.23，趋于 ∞ 时 1.96）：自由度越大 → 分母噪声越小 → 尾部越薄 → 临界值越小。详见 [[hypothesis_testing]]。

## Related Notes

- [[hypothesis_testing]]
