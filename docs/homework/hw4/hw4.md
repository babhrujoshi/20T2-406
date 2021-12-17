# **CPSC 406: Homework 4 (Due Mar 26, 11:59pm)**


1. Consider 

    \begin{equation}\label{e-constrained}
    \begin{array}{ll}
    \min_{x}& f(x)\\
    \st & x\in \mathcal C
    \end{array}
    \end{equation}

    where $f(x)$ is a convex and differentiable function, and $\mathcal C$ is a convex set. We know that the projected gradient method converges to an optimal solution, so we know that $x^*$ is optimal if and only if 

    \begin{equation}\label{e-projection}
    x^* = \proj_{\mathcal C}(x^* - \gamma \nabla f(x^*))
    \end{equation}

    for any constant $\gamma > 0$, where the projection itself is an optimization problem:

    $$
    \begin{array}{lll}
    \proj_{\mathcal C}(z) = &\underset{x}{\text{argmin}}
    & \|x - z\|_2^2\\
    &\st & x\in \mathcal C
    \end{array}
    $$

    1. **Nonnegative constraint.** For $\mathcal C = \{x : x_i \geq 0\}$,  prove that \eqref{e-projection} implies $\nabla f(x^*)\geq 0$.
    
    2. **Normal cone.** In general, the optimality condition for \eqref{e-constrained} is that $x^*$ is optimal if and only if
    
    \begin{equation}\label{e-optimality}
    \nabla f(x^*)^T(y-x^*) \geq 0,\; \forall y\in \mathcal C
    \end{equation}

    In other words, $-\nabla f(x^*)$ is in the *normal cone of $\mathcal C$ at point $x^*$*. Show that this property is equivalent to \eqref{e-projection}.

  

2.  Let $\mathcal{C} = \{x \in \R^n | Ax = b\}$, where $A \in \R^{m\times n}$ has full row rank. Show that $\proj_{\mathcal C} (z) = z + A\trans(AA\trans)^{-1}(b-Az)$. Hint: decompose $\nabla f(x)$ into two components: $\nabla f(x)$ is in $\range(A^T) \oplus \vnull(A)$.

3. **Projected gradient descent** Consider the minimization problem

    \begin{align*}
    &\min 2x_1^2 + 3x_2^2 + 4x_3^2 + 2 x_1 x_2 -2x_1 x_3 - 8x_1 - 4x_2 - 2x_3\\
    & \text{subject to } x_1, x_2, x_3 \geq 0
    \end{align*}

    1. Show that the vector $(\frac{17}{7},0,\frac{6}{7})^\intercal$ is an optimal solution of above minimization problem.

    2. Employ gradient projection method with constant stepsize $\frac{1}{L}$ ($L$ being the Lipschitz constant of the gradient of the objective function). Show the function value and the solution after the first 100 iteration.

    3. Using [Convex.jl](https://jump.dev/Convex.jl/v0.13.2/index.html) or [CVX](http://cvxr.com/cvx/) package solve the above minimization program and compare the solution with the solution from part b.

4. **Convexity**

    1. **Entropy** Let $f(x) = -\sum_{i=1}^n x_i\log x_i$. Show that $f(x)$ is concave in the simplex set $\{ x\in \R^n | \sum_{i=1}^n x_i = 1, x_1,\dots,x_n\geq 0 \}$.

    2. **Log-Sum-Exp** Let $A \in \R^{m\times n}$ with $a_i$ as the $i$th column of $A\trans$. Show that $f(x) = \log\left(\sum_{i=1}^m e^{a_i\trans x}\right)$ is convex in $\R^n$.


