# **CPSC 406: Homework 3 (Due March 10, 11:59pm)**


1.  Consider $\mA\in \R^{m\times n}$ and $\mB\in \R^{m\times n}$. Show that 
    $$ \tr(\mA\trans\mB) = \sum_{i,j} A_{ij}B_{ij}.$$ 
    Compute the number of flops and storage required to calculate the left hand side ($\tr(\mA'*\mB)$) and right hand side ($\sum_{i,j} A_{ij}B_{ij}$) . (Use  $O(f(n,m))$ notation.) Which of the two is the "superior" formulation?

2.  (Structured matrices.) Consider a  matrix $\mA\in \R^{n\times n}$.

    1.  Suppose $\mA$ is diagonal, e.g. $A_{ij} = 0$ whenever $i\neq j$. Prove that $\mA \succeq 0$ if and only if $A_{ii} \geq 0$ for all $i$, and $\mA\succ 0$ if and only if $A_{ii} > 0$ for all $i$.

    2.  Suppose that $\mA$ is block diagonal, e.g. 

        $$ \mA = \bmat \mA^{(1)} & 0 &  \cdots & 0 \\ 0 & \mA^{(2)} & \cdots & 0 \\ \vdots   & \vdots   & \ddots & \vdots \\0 & 0 & \cdots  &\mA^{(l)}  \\ \emat.$$

        Prove that $\mA\succeq 0$ if and only if each block $\mA^{(k)} \succeq 0$ for all $k$, and $\mA\succ 0$ if and only if $\mA^{(k)} \succ 0$ for all $k$.

3.  We will show that if a twice continuously differentiable $f:\R^n\rightarrow \R$ has $L$-Lipschitz continous gradient and a minimizer of $\min_{\vx \in \R^n} f(\vx)$ exists, then the graident method with constant step size converges if step size $\bar{\alpha} \in (0,\frac{2}{L})$. In the following subquestion, assume that $f$ has an $L$-Lipschitz continous gradient and is twice continuously differentiable.

    1.  Fix $\vx \in \R^n$. Show that $\|\nabla^2 f(\vx)\vv \|_2 \leq L\|\vv\|_2$ for all $\vv$. (Hint: Let $c:\R^n\rightarrow\R^n$ be continuously differentiable. Then the directional derivative of $c$ at $\vx$ in the direction of $\vv$ is 

        $$ \mJ\vv = \lim_{t ↘ 0} \frac{c(\vx + t\vv) - c(\vx)}{t},$$

        where $\mJ$ is the Jacobian of $c$ at $\vx$. Apply this to $c = \nabla f$ and take 2-norm of both sides.)

    2.  Show that for all $\vx$, the eigenvalues of $\nabla^2 f(\vx)$ are bounded from above by $L$ , i.e. for all $\vx$, we have $\nabla^2 f(\vx) \preceq L \mI$. 

    3.  Building on part b. and using multivariate Taylor's remainder theorem (Beck Thm 1.24), show that

        $$f(\vv) \leq f(\vw) + \nabla f(\vw)\trans (\vv-\vw) + \frac{L}{2}\|\vv-\vw\|_2^2$$

        for all $\vv$ and $\vw$. This is known as the descent lemma.

    4.  Consider gradient descent with constant step size $\alpha$, i.e.

        $$\vx_{k+1} = \vx_k - \alpha \nabla f(\vx_k).$$

        In the descent lemma, substitute $\vv = \vx_{k+1}$ and $\vw = \vx_{k}$ and show that $f(\vx_{k+1}) < f(\vx_k)$ if $\alpha \in (0,\frac{2}{L})$ and $\vx_k$ is not a stationary point.

5. Prove that the intersection of convex sets 

    $$
    \mathcal{S} = \mathcal{S}_1\cap \mathcal{S}_2 \cap \cdots \cap \mathcal{S}_n = \{x : x\in \mathcal{S}_1 \text{ and } \mathcal{S}_2 \text{ and } \cdots \text{ and } \mathcal{S}_n\}
    $$

    is itself a convex set.

6. The *second order cone* is the set of tuples $(x\in \R^n, t\in \R)$ 

    $$
    \mathcal{S} = \{(x,t) : \|x\|_2 \leq t\}.
    $$

    Prove that $\mathcal{S}$ is convex.