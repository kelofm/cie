# PyFEM

This module is meant to be a framework for supporting proofs-of-concept of the FEM-based optimal control of PDEs (especially the heat equation) using adjoint formulations. To put it mildly, performance is not a primary aspect of this framework.


# FE Formulation

## Stationary Heat Conduction

Let's start with a 1D stationary nonlinear heat conduction with Neumann boundaries:
$$
\frac{d}{dx}( \kappa \frac{du}{dx} ) = \lambda f                \\
\frac{du}{dx}\bigg\rvert_{x_0} = \lambda \Phi_0                 \\
\frac{du}{dx}\bigg\rvert_{x_1} = \lambda \Phi_1
\tag{1}
 $$
* $`u = u(x)`$ temperature field
* $`\kappa = \kappa (u(x))`$ conductivity
* $`f=f(u(x))`$ load
* $`\lambda`$ load factor

Unpacking the derivative using the chain rule yields:
$$
\frac{d\kappa}{dx}\frac{du}{dx} + \kappa \frac{d^2u}{dx^2}
- \lambda f
= 0         \tag{2}
$$


$$
\frac{d\kappa}{du} (\frac{du}{dx})^2 + \kappa \frac{d^2u}{dx^2}
- \lambda f
= 0         \tag{3}
 $$

Choose an appropriate function space $`V`$ and assume that the solution $`u = \sum\limits_i N_i \hat{u}_i \in V`$ and the test function $`v = \sum\limits_i N_i \in V`$, then multiply by the test function $`v`$ and integrate over the domain $`\Omega=[x_0, x_1]`$:

$$
\int_\Omega v \frac{d\kappa}{du} (\frac{du}{dx})^2 d\Omega + \int_\Omega \kappa v \frac{d^2u}{dx^2} d\Omega
- \int_\Omega \lambda v f d\Omega
= 0         \tag{4}
 $$

To sneak in the boundaries and get rid of the second derivatives, the second term on the LHS of $`(5)`$ has to be partially integrated (a threeway partial integration):
$$
\int_\Omega v \frac{d\kappa}{du} (\frac{du}{dx})^2 d\Omega
+ \bigg[\kappa v \frac{du}{dx} \bigg]_{x_0}^{x_1}
- \int_\Omega v \frac{d\kappa}{du} (\frac{du}{dx})^2 d\Omega
- \int_\Omega \kappa \frac{dv}{dx} \frac{du}{dx}d\Omega
- \lambda \int_\Omega v f d\Omega
= 0         \tag{5}
 $$

$$
\int_\Omega \kappa \frac{dv}{dx} \frac{du}{dx}d\Omega
- \lambda (\kappa_1 v_1 \Phi_1 - \kappa_0 v_0 \Phi_0)
+ \lambda \int_\Omega v f d\Omega
= 0         \tag{6}
 $$

* $`\kappa_0=\kappa(u(x_0))`$
* $`\kappa_1=\kappa(u(x_1))`$

The solutions at the boundaries $`u_0, u_1`$ are of course still unknown. Substituting the discrete counterparts of $`u`$ and $`v`$  in $`(6)`$ yields:

$$
\int_\Omega \kappa \frac{dN_i}{dx} \frac{dN_j}{dx} d\Omega \hat{u}_i
- \lambda (\kappa_1 N_i(x_1) \Phi_1 - \kappa_0 N_i(x_0) \Phi_0)
+ \lambda \int_\Omega f N_i
= 0         \tag{7}
 $$

At this point, the stiffness matrix and the load vector could be computed if $`\kappa`$ was constant, and the system could be solved. However, we're focusing on the nonlinear system and have to proceed with some nonlinear theory.

From now on, the LHS in $`(7)`$ will be referred to as the residual $`r=r(u,\lambda)`$, and will be the basis of the analysis. The problem boils down to finding the roots of this residual for specified load factors, which is done by a combination of explicit ODE integration (predictor) and some variation of the Newton method (corrector). Unfortunately, we'll need the derivatives for the Newton method.

The partial derivatives of the residual form the tangent stiffness matrix $`\bar K`$ and load vector $`q`$:

$$
\bar{K}_{ij} = \frac{\partial r_i}{\partial \hat{u}_j}
=
\int_\Omega \kappa \frac{dN_i}{dx}\frac{dN_j}{dx}d\Omega
+ \int_\Omega \frac{d\kappa}{du}N_j\frac{dN_i}{dx}\frac{dN_k}{dx}d\Omega \hat{u}_k
+ \lambda \int_\Omega \frac{\partial f}{\partial \hat u_j}N_i N_j d\Omega
\tag{8}
 $$

$$
q_i = -\frac{\partial r_i}{\partial \lambda}
=
- \int_\Omega f N_i d\Omega + \bigg( \kappa_1 N_i(x_1) \Phi_1 - \kappa_0 N_i(x_0) \Phi_0 \bigg)
\tag{9}
 $$

For the sake if simplicity, the load $`f`$ will be assumed to be independent of $`u`$. Then, the **tangent stiffness matrix $`\bar{K}`$** can be broken down into two components: the **stiffness matrix $`K`$** and the **geometric stiffness matrix $`K_g`$**.

$$
\bar{K} = K + K_g
\tag{10}
 $$

* $`K       = \int_\Omega \kappa \frac{dN_i}{dx}\frac{dN_j}{dx}d\Omega`$
* $`K_g     = \int_\Omega \frac{d\kappa}{du}N_j\frac{dN_i}{dx}\frac{dN_k}{dx} d\Omega \hat{u}_k  `$

## Alternative to the Newton method

Since the Newton method comes with a lot of extra computational burdens while not ensuring convergence, it's a good idea to try cheaper alternatives. One possible method is converting the problem into a fixed point iteration:

$$
\Phi(\hat u)
=
K(\hat u)\hat u - f(\hat u)
= 0     \tag{11}
 $$

Of course it might not converge either, but at least we don't have to compute the geometric stiffness $`K_g`$ for example. Also, the iteration might converge to a different solution, but we can detect most of those cases by keeping an eye on the residual, which is supposed to approach $`0`$.

The most scary situation would be converging to a solution that **is** in equilibrium ($`r=0`$), but on a different equilibrium path, which would indicate bifurcations or very tricky nonlinearities. Thankfully, we don't have to worry about bifurcations for heat problems (to the extent of my knowledge).


## Transient Heat Conduction

The formulation of the transient 1D nonlinear heat conduction will be mostly identical to the stationary case, but with the added time derivative term:

$$
-c \frac{\partial u}{\partial t} + \frac{\partial }{\partial x}(\kappa \frac{\partial u}{\partial x}) - \lambda f
= 0
\\
\frac{\partial u}{\partial x} \bigg \rvert_{x_0} = \lambda \Phi_0(t)
\\
\frac{\partial u}{\partial x} \bigg \rvert_{x_1} = \lambda \Phi_1(t)
\\
u(x,0) = u_0(x)
\tag{12}
 $$

* $`u = u(x, t)`$
* $`c = c(u)`$
* $`\kappa = \kappa(u)`$

To reduce clutter, the following notation is intruduced:
* $`\dot \square    := \frac{\partial \square}{\partial t}`$
* $`\square'        := \frac{\partial \square}{\partial x}`$

Expand the derivatives:
$$
-c \dot u + \frac{d\kappa}{du} (u')^2 + \kappa u'' - \lambda f
= 0
\tag{13}
 $$

Multiply with the test function $`v`$ and integrate over the space domain:
$$
-\int_\Omega c v \dot u d\Omega
+ \int_\Omega \frac{d\kappa}{du} v (u')^2 d\Omega
+ \int_\Omega \kappa v u'' d\Omega
- \lambda \int_\Omega f v d\Omega
= 0
\tag{14}
 $$

Integrate by parts:
$$
-\int_\Omega c v \dot u d\Omega
+ \bigg[ \kappa v u' \bigg]_{x_0}^{x_1}
- \int_\Omega \kappa v' u' d\Omega
- \lambda \int_\Omega f v d\Omega
= 0
\tag{15}
 $$

Substituting the discrete counterparts of $`u`$ and $`v`$ works the same way as in the stationary case, except the coefficients are now time-dependent $`\hat u_i = \hat u_i (t)`$.
$$
\int_\Omega c N_i N_j d\Omega \dot{\hat u}_i
+ \int_\Omega \kappa N_i' N_j' d\Omega \hat u_i
- \lambda \bigg( \kappa_1N_i(x_1)\Phi_1 - \kappa_0 N_i(x_0) \Phi_0 \bigg)
+ \lambda \int_\Omega f N_i d\Omega
= 0
\tag{16}
 $$

Before formulating the residual $`r(u,\lambda)`$, we have to discretize in time as well. However, the current notation is a bit hard to work with so from this point, I'm just going to use the universal notation for "structural" matrices and the load vector:
* $`M   := \int_\Omega c N_i N_j d\Omega`$
* $`K   := \int_\Omega \kappa N_i' N_j' d\Omega`$
* $`q   := - \bigg( \kappa_1N_i(x_1)\Phi_1 - \kappa_0 N_i(x_0) \Phi_0 \bigg) + \int_\Omega f N_i d\Omega`$

Keep in mind that these entities are not truly derived from the residuals yet. For example, the "stiffness matrix" $`K_{ij}`$ is obviously not equal to $`\frac{\partial r_i}{\partial \hat u_j}`$ due to the non-constant material properties $`c(u)`$ and $`\kappa(u)`$, but we'll use them for the time being to keep things readable.

Now we can discretize in time as well, which we'll do using a mixed explicit-implicit finite difference scheme. The basic idea is as follows:
$$
\frac{dF(t)}{dt} = G(t,F(t))
\rightarrow
\frac{F^{k+1} - F^{k} }{ \Delta t }
\approx
\theta G(t^{k+1}, F^{k+1}) + (1-\theta) G(t^k, F^k)
\tag{17}
 $$

* $`\Delta t := t^{k+1} - t^k`$
* $`\square^k := \square(t^k)`$
* $`\theta \in [0,1]`$


Apply $`(17)`$ to $`(16)`$:
$$
\frac{\hat u^{k+1} - \hat u^k}{\Delta t}
+ \theta (M^{k+1})^{^{-1}} \bigg( K^{k+1} \hat u^{k+1} - \lambda^{k+1} q^{k+1} \bigg)
+ (1-\theta)M^{k^{-1}} \bigg( K^k \hat u^k - \lambda^k q^k  \bigg)
= 0
\tag{18}
 $$

To sort out this mess, we need to separate the unknown terms ($`\square^{k+1}`$) from the terms that we already computed ($`\square^k`$):
$$
\bigg[
    \bigg( \theta (M^{k+1})^{^{-1}} K^{k+1} + \frac{1}{\Delta t} I \bigg) \hat u^{k+1} - \theta (M^{k+1})^{^{-1}} \bigg( \lambda^{k+1} q^{k+1} \bigg)
\bigg]
+
\bigg[
    (1-\theta) M^{k^{-1}} \bigg( K^k \hat u^k - \lambda^k q^k \bigg)
    - \frac{1}{\Delta t} \hat u^k
\bigg]
= 0
\tag{19}
 $$

*Several things to note here:*
*First of all, it might seem weird to carry the load factors ($`\lambda`$) through the time steps but if we don't want to impose extra constraints on the nonlinear iteration later, we have to keep using them. The other painful thing is that the mass matrix ($`M`$) has to be inverted at every time step and every nonlinear iteration step. This is a huge difference compared to the solution of the linear system, where the mass matrix is constant, and can just premultiply the affected terms without ever having to invert it. As reference, here's how $`(19)`$ looks like in the linear case (premultiplied by $`M`$):*
$$
\bigg[
    \bigg( \theta K + \frac{1}{\Delta t}M \bigg) \hat u^{k+1}
\bigg]
+
\bigg[
    \bigg( (1-\theta) K - \frac{1}{\Delta t} M \bigg)
    \hat u^k
    - (1-\theta) q^k
    - \theta q^{k+1}
\bigg]
= 0
\tag{20}
 $$


Back from sidetracking, the task now is to compute the derivative of the transient residual $`r^{k+1}`$ $`(19)`$ with respect to $`\hat u^{k+1}`$. Avoiding an infinite slew of indices, the following expressions are written in another kind of notational catastrophe: mixed matrix and index notation.
$$
\frac{\partial r^{k+1}}{\partial \hat u^{k+1}_l}
=
\frac{\partial}{\partial \hat u^{k+1}_l}
\bigg(
    \theta M^{k+1^{-1}} K^{k+1} \hat u^{k+1}
    +
    \frac{\hat u^{k+1}}{\Delta t}
    -
    \theta M^{k+1^{-1}} \lambda^{k+1} q^{k+1}
\bigg)
\tag{21}
 $$

From this point the time index is dropped since all components live in the $`(k+1)^{th}`$ time step. Rearranging the expression above and highlighting important terms:
$$
\frac{\partial r}{\partial \hat u_l}
=
\theta \textcolor{red}{\frac{\partial M^{-1}}{\partial \hat u_l} K \hat u}
+
\theta M^{-1} \textcolor{green}{\frac{\partial (K \hat u) }{\partial \hat u_l}}
+
\frac{1}{\Delta t}\frac{\partial \hat u}{\partial \hat u_l}
-
\lambda \theta \textcolor{teal}{\frac{\partial M^{-1}}{\partial \hat u_l} q}
\tag{22}
 $$

The expression in $`\textcolor{green}{green}`$ should be familiar; it's the tangent stiffness matrix $`\bar K`$ from the stationary case $`(10)`$. The other two coloured expressions are similar, in that they are products of $`3^{rd}`$ order tensors and vectors, resulting in $`2^{nd}`$ order tensors. Expanding the interesting bits of $`(22)`$, we get:
$$
\frac{\partial r_i}{\partial \hat u_l}
=
-\theta \textcolor{red}{
    [M^{-1}]_{ip} \int_\Omega \frac{dc}{du} N_p N_k N_j d\Omega [M^{-1} K \hat u]_j
    }
+
\theta [M^{-1} (\textcolor{green}{K+K_g}) ]_{il}
+
\frac{1}{\Delta t} \delta_{il}
+
\lambda \theta \textcolor{teal}{
    [M^{-1}]_{ip} \int_\Omega \frac{dc}{du} N_p N_k N_j d\Omega [M^{-1} q]_j
    }
\tag{23}
 $$

And we're done. The residual $`(19)`$ and its derivative $`(22)`$ can be used for iteratively solving the nonlinear equation for each successive time step.


## Adjoint

Based on Julian's thesis, the adjoint of the nonlinear transient heat equation takes the following form:
$$
-
\bigg(
    \frac{\partial c}{\partial u} \bigg \vert_{\bar u} \bar u
    + c(\bar u)
\bigg) \partial_\tau p
+
\bigg(
    \frac{ \partial \kappa }{ \partial u} \bigg \vert_{\bar u }
    + \kappa (\bar u)
\bigg) \Delta p
=
u_\Omega -\bar u
\tag{24}
 $$


$$
\bigg(
    \kappa(\bar u)
    + \frac{ \partial \kappa }{ \partial u } \bigg \vert_{\bar u} \bar u
\bigg) \partial_n p
+
\bigg(
    \frac{ \partial \kappa } { \partial u } \bigg \vert_{\bar u} \bar u
    - \alpha \kappa(\bar u) p
\bigg)
=
0
\tag{25}
 $$


$$
p(\bullet, 0) = 0
\tag{26}
 $$

* $`p = p(\tau,x)`$ : adjoint state
* $`\tau = T - t`$ : reversed time variable
* $`\bar u`$ : forward solution
* $`u_\Omega`$ : prescribed state

*Note:*
*The initial condition $`(26)`$ is most likely incorrect, as it would mean that no control is needed at $`\tau=0`$ ($`t=T`$)*.

Keep in mind that the time variable has already been reversed here.

The DE $`(24)`$ is almost identical to the original one $`(12)`$, save for the coefficients of the adjoint state $`p`$, and the space derivatives. The coefficients are now functions of $`\tau`$ and $`x`$ but not $`p`$, which makes the derivation a bit simpler.

Let's rename the coefficients from $`(24)`$ and reduce the system to 1D:
$$
- \tilde c \dot p + \tilde \kappa p'' = \tilde f
\tag{27}
 $$

$$
\tilde c(x) := \frac{\partial c}{\partial u} \bigg \vert_{\bar u} \bar u + c(\bar u)
\tag{28}
 $$

$$
\tilde \kappa(x) := \frac{ \partial \kappa }{ \partial u} \bigg \vert_{\bar u } \bar u + \kappa (\bar u)
\tag{29}
 $$

$$
\tilde f(x) := u_\Omega - \bar u
\tag{30}
 $$

Note that the time derivative is replaced with the derivative of the reversed time variable $`\dot p := \frac{\partial p}{\partial \tau}`$

Following the same line of thought as before:
$$
- \int_\Omega \tilde c v \dot p d \Omega
+
\int_\Omega \tilde \kappa v p'' d\Omega
=
\int_\Omega v \tilde f d\Omega
\tag{31}
 $$

However, integrating by parts yields an extra (non-symmetric) term because $`\tilde \kappa`$ does not depend on $`p`$ :
$$
-\int_\Omega \tilde c v \dot p d\Omega
+
\bigg[ \tilde \kappa v p' \bigg]_{\partial \Omega}
-
\int_\Omega \tilde \kappa' v p' d\Omega
-
\int_\Omega \tilde \kappa v' p' d\Omega
-
\int v \tilde f d\Omega
=
0
\tag{32}
 $$

$$
-\int_\Omega \tilde c N_i N_j d\Omega \dot{\hat p}_j
+
\bigg[ \tilde \kappa v p' \bigg]_{\partial \Omega}
-
\int_\Omega \tilde \kappa' N_i N_j'd\Omega \hat p_j
-
\int_\Omega \tilde \kappa N_i' N_j' d\Omega \hat p_j
-
\int_\Omega N_i \tilde f d\Omega
=
0
\tag{33}
 $$

Similarly introducing the *mass* and *stiffness* matrices, the boundary term $`b_i`$, the load $`\mathfrak{f}_i`$ , plus the new matrix $`L_{ij}`$ resulting from the non-symmetric term:
$$
-M_{ij}\dot{\hat p}_j + b_i - \bigg( L_{ij} + K_{ij} \bigg) \hat p_j - \mathfrak{f}_i
=
0
\tag{34}
 $$


Now we can discretize in time:
$$
\dot{\hat p}^{k+1} \approx \frac{ \hat p^{k+1} - \hat p^k }{\Delta t}
=
\theta M^{k+1^{-1}}
\bigg[
    b^{k+1} - (K^{k+1} + L^{k+1}) \hat p^{k+1} - f^{k+1}
\bigg]
+
(1-\theta) M^{k^{-1}}
\bigg[
    b^k - (K^k + L^k) p^k - f^k
\bigg]
\tag{35}
 $$

$$
\bigg[
    \frac{1}{\Delta t} I + \theta M^{k+1^{-1}} (K^{k+1} + L^{k+1})
\bigg] \hat p^{k+1}
=
\frac{1}{\Delta t} \hat p^k
+
\theta M^{k+1^{-1}} (b^{k+1}-f^{k+1})
+
(1-\theta) M^{k^{-1}}
\bigg[
    b^k - (K^k + L^k) \hat p^k - f^k
\bigg]
\tag{36}
 $$

Nothing depends on $`\hat p^{k+1}`$, so we don't need the derivative of $`(36)`$ (which is why no load factor $`\lambda`$ was included in the derivation).


# Module

## Class Structure

In most cases, I tried to create an object-oriented structure but you might find some exceptions. Anyway, here's brief overview:

* [*Integrator*](https://gitlab.com/matekelemen/CiE/-/tree/optcontrol/libraries/fem/python/modules/pyfem/numeric/integration.py) is responsible for integrating a scalar, single variate input function on a specified domain. The current implementation uses gaussian quadrature with user-selectable polynomial order.

* [*BasisFunctions*](https://gitlab.com/matekelemen/CiE/-/tree/optcontrol/libraries/fem/python/modules/pyfem/discretization/basisfunctions.py): 1D matching set of basis functions defined on the [-1,1] domain. The functions in this class form a basis for their respective function spaces. *BasisFunctions* is responsible for providing an interface for the evaluation of its functions at one or several points, and evaluating to 0 when sampled out of bounds. Currently, all scripts use integrated legendre basis functions (even though a polynomial order of 1 is used, meaning that these basis functions are identical to the classic linear ones).

* [*Element*](https://gitlab.com/matekelemen/CiE/-/tree/optcontrol/libraries/fem/python/modules/pyfem/discretization/element.py) has instances of *BasisFunctions* (1 in 1D) that it maps to its own domain, a list of degrees-of-freedom that maps their local DoFs to the global DoFs the *Element* belongs to, and an instance of *Integrator*. *Element* is responsible for computing its own stiffness matrix, mass matrix, and load vector, as well as updating the global ones. Also, it is responsible for computing the solution field (temperature field for example) in its own domain by evaluating its basis functions weighted with coefficients.

* [*BoundaryCondition*](https://gitlab.com/matekelemen/CiE/-/tree/optcontrol/libraries/fem/python/modules/pyfem/discretization/boundaryconditions.py) might not be the best abstraction, so it may be removed in the future. Currently its only function is to store the global DoF it acts on, and the value of the boundary condition (and also a penalty value depending on the type of the boundary condition). Since different physical models require different enforcement of boundary conditions, this task is delegated to *FEModel*.

* [*FEModel*](https://gitlab.com/matekelemen/CiE/-/tree/optcontrol/libraries/fem/python/modules/pyfem/discretization/femodel.py) has a list of *Element*s and *BoundaryCondition*s, and is responsible for computing the global structural matrices and the load vector, along with modifying them to enforce boundary conditions (*TransientFEModel*s have time-dependent boundaries and time-and-space-dependent loads). *FEModel* is also responsible for evaluating the solution field (heat field for example), given a set of coefficients for the global DoFs, by evaluating the individual contributions of all elements on their respective domains. Note that *FEModel* **is not responsible for the assembly** process, meaning that you have to come up with the distribution of global DoFs yourself and pass subsets of it to *Element*s such that they are connected the way you want them to be.


This concludes the basic structure of classes for the FE framework, which might be expanded in the future.