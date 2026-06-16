

## 1 策略是否失效问题
Probability an asset with annual Sharpe ratio 1 loses money in 4 years.

The annual Sharpe ratio is given as $1$, meaning $\frac{\mu}{\sigma} = 1$, which simplifies to **$\mu = \sigma$**.

$$\text{Mean}_{4\text{yr}} = 4\mu$$$$\text{Volatility}_{4\text{yr}} = \sqrt{4\sigma^2} = 2\sigma$$We want to find the probability that the total 4-year return is less than zero. We can do this by calculating the Z-score, which measures how many standard deviations zero is away from the expected 4-year return.

$$Z = \frac{\text{Target} - \text{Mean}_{4\text{yr}}}{\text{Volatility}_{4\text{yr}}}$$

$$Z = \frac{0 - 4\mu}{2\sigma}$$$$Z = \frac{-4\mu}{2\mu} = -2$$According to standard normal distribution tables, **Φ(−2)≈0.025**, meaning there is roughly a **2.5%** chance the asset loses money over a 4-year stretch.


## 2 MM豆的统计检验
A person claims they can identify the color of M&M candies by taste. There are 5
colors. In a test of 3 candies, they get 2 correct. What can you conclude about the claim？

To prove this, we calculate the p-value (the probability of getting 2 or more correct by random guessing) using a Binomial distribution where $n=3$ and the probability of a correct guess is $p=0.2$.

- **Probability of exactly 2 correct:**
    
    $$P(X=2) = \binom{3}{2} (0.2)^2 (0.8)^1 = 3 \times 0.04 \times 0.8 = 0.096$$
    
- **Probability of exactly 3 correct:**
    
    $$P(X=3) = \binom{3}{3} (0.2)^3 (0.8)^0 = 1 \times 0.008 \times 1 = 0.008$$
    

The total probability of getting at least 2 correct is $0.096 + 0.008 = \mathbf{0.104}$.

Because there is a **$10.4\%$ chance** of this happening by pure luck, it exceeds the standard statistical significance threshold ($\alpha = 0.05$). We fail to reject the null hypothesis; there is not enough evidence to support their claim.

(a) Follow-up: What if the person tastes 100 candies and gets 40 correct? (Use an
appropriate approximation.)

clearly, this is sum of 100 Bernouli experiment with p =1/5.
$SE(p) = \sqrt{\frac{p(1-p)}{n}}$，**在 M&M 例子中：** $n=100$, $p=0.2$。所以 $SE = \sqrt{\frac{0.2 \times 0.8}{100}} = \sqrt{0.0016} = 0.04$ （也就是 **4%**. CI [0.12,0.28], 0.4>0.28 -> I believe he has the ability

注意，在构建置信区间的时候，显而易见有两种可能性，比如上面的问题，p可以用原假设的0.2，也可以用通过观测值推算的 $\hat p = 40/100 = 0.4$, 在做假设检验的时候，我们一般以原假设的p构建接受区间（acceptance region）,如果观测点落在acceptance region 外 (一般是 正负2个标准差$se(\hat p))$，我们拒绝原假设。而如果用p_hat构建区间，照理说也是可以的，如果原假设落在p_hat构建的区间外，我们应该也可以拒绝原假设，但保险来说，我们如果有原假设的参数量，尽量用原假设下的参数构建Z 统计量，或者T统计量，或者原假设下的接受区间。

(b) Follow-up: If they are not randomly guessing yet average 4 correct per 10, provide a
reasonable explanation for this pattern.

**Specific Dye Sensitivity:** They might only be able to taste one or two specific colors. For example, if red and green dyes have a distinct bitter chemical taste to them, they can identify those with 100% accuracy, but they must guess randomly among the remaining three colors. This combination of perfect knowledge for some and guessing for others would elevate their overall average above 20%.

(c) Follow-up: The choice of $\alpha$ and Mind Reading

For a claim like mind reading, you should choose an incredibly strict, microscopic $\alpha$ (e.g., $\alpha = 0.0000003$, which corresponds to a $5\sigma$ standard used in particle physics).
This is rooted in Bayesian statistics and the famous adage: _"Extraordinary claims require extraordinary evidence."_ * A claim of mind reading violates all known laws of physics. Therefore, the **prior probability** of it being true is effectively zero.