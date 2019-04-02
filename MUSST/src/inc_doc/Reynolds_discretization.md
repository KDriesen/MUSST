**REYNOLDS EQUATION DISCRETIZATION**

<!-- ========== Reynolds equation discretization ========== -->

<div><button class="collapsible">Reynolds equation discretization</button>
<div class="content">

\(
   \def\scall{{ \, \tiny{\bullet} \, }}
\)

The Reynolds' equation writes:

\[
\text{div}\left( \frac{\rho h^3}{\mu} \overrightarrow{\text{grad}} \, p \right) = 6 \, \text{div}\left( \rho h \, \vec{u} \right)
\]

Over each element \( \Omega_e \), the equation is multiplied by the weighting function \( W_i \).

Owing the relationship:

\[
b \, \text{div}\vec{a} = \text{div} (b\,\vec{a}) - \vec{a} \scall \overrightarrow{\text{grad}}b
\]

the equation becomes:

\[
\begin{split}
\int \!\!\! \int_{\Omega_e} \text{div} \left( \frac{\rho h^3}{\mu}  W_i \, \overrightarrow{\text{grad}} \, p \right) d\Omega_e &-
\int \!\!\! \int_{\Omega_e} \frac{\rho h^3}{\mu} \overrightarrow{\text{grad}} \, p \, \scall \, \overrightarrow{\text{grad}} \, W_i \, \mathrm{d}\Omega_e = \\
\int \!\!\! \int_{\Omega_e} 6 \, \text{div} \left( \rho h  W_i \, \vec{u} \right) \, \mathrm{d}\Omega_e &-
\int \!\!\! \int_{\Omega_e} 6 \, \rho h  \, \vec{u} \, \scall \, \overrightarrow{\text{grad}} \, W_i \, \mathrm{d}\Omega_e
\end{split}
\]

\( W_i \) vanishes on \( \Omega_e \) frontier \( \Gamma_e \), then according Green's formula:

\[
\int \!\!\! \int_{\Omega_e} \text{div} \left( \frac{\rho h^3}{\mu}  W_i \, \overrightarrow{\text{grad}} \, p \right) \mathrm{d}\Omega_e =
\oint_{\Gamma_e} \frac{\rho h^3}{\mu}  W_i \, \overrightarrow{\text{grad}} \, p \scall \vec{n} \, \mathrm{d}\Gamma_e = 0
\]

In the same manner, if the right handside is handled likewise:

\[
\int \!\!\! \int_{\Omega_e} 6 \, \text{div} \left( \rho h  W_i \, \vec{u} \right) \, \mathrm{d}\Omega_e =
\oint_{\Gamma_e} 6 \, \rho h  W_i \, \vec{u} \scall \vec{n} \, \mathrm{d}\Gamma_e = 0
\]

and then:

\begin{equation} \label{trans}
\int \!\!\! \int_{\Omega_e} \frac{\rho h^3}{\mu} \overrightarrow{\text{grad}} \, p \scall \overrightarrow{\text{grad}} \, W_i \, \mathrm{d}\Omega_e =
\int \!\!\! \int_{\Omega_e} 6 \, \rho h  \, \vec{u} \scall \overrightarrow{\text{grad}} \, W_i \, \mathrm{d}\Omega_e
\end{equation}

If the right handside isn't transformed:

\begin{equation} \label{notrans}
\int \!\!\! \int_{\Omega_e} \frac{\rho h^3}{\mu} \overrightarrow{\text{grad}} \, p \scall \overrightarrow{\text{grad}} \, W_i \, \mathrm{d}\Omega_e =
\int \!\!\! \int_{\Omega_e} 6 \, \text{div}\left( \rho h \, \vec{u} \right) \, W_i \, \mathrm{d}\Omega_e
\end{equation}

Equation \eqref{notrans} is more suitable when linear shape functions are use with SUPG, otherwise, using equation \eqref{trans}, \( \overrightarrow{\text{grad}} \, W_i = \overrightarrow{\text{grad}} \, N_i \).

As:

\begin{align*}
\int \!\!\! \int_{\Omega_e} 6 \, \text{div}\left( \rho h \, \vec{u} \right) \, N_i \, \mathrm{d}\Omega_e & = -\int \!\!\! \int_{\Omega_e} 6 \, \rho h  \, \vec{u} \scall \overrightarrow{\text{grad}} \, N_i \, \mathrm{d}\Omega_e \\
\vec{\alpha} & = a \, \vec{u}
\end{align*}

then:

\begin{align*}
\int \!\!\! \int_{\Omega_e} 6 \, \text{div}\left( \rho h \, \vec{u} \right) \, W_i \, \mathrm{d}\Omega_e
& = + \int \!\!\! \int_{\Omega_e} 6 \, \text{div}\left( \rho h \, \vec{u} \right) \, \left( N_i + \vec{\alpha} \scall \overrightarrow{\text{grad}} \, N_i \right) \mathrm{d}\Omega_e \\
& = - \int \!\!\! \int_{\Omega_e} 6 \, \rho h  \, \vec{u} \scall \overrightarrow{\text{grad}} \, N_i \, \mathrm{d}\Omega_e \\
&  \phantom{ = } + \int \!\!\! \int_{\Omega_e} 6 \, \text{div}\left( \rho h \, \vec{u} \right) \, \vec{\alpha} \scall \overrightarrow{\text{grad}} \, N_i \, \mathrm{d}\Omega_e \\
& = -\int \!\!\! \int_{\Omega_e} 6 \, \vec{u} \scall \overrightarrow{\text{grad}} \, N_i \left( \rho h - \vec{\alpha} \scall \overrightarrow{\text{grad}} \, \rho h \right) \mathrm{d}\Omega_e
\end{align*}

Approximating \( \rho h \) with the shape functions, \( \rho h = [\rho h]_k \, N_k \), leads to:

\begin{align*}
\int \!\!\! \int_{\Omega_e} 6 \, \text{div}\left( \rho h \, \vec{u} \right) \, W_i \, \mathrm{d}\Omega_e & = -\int \!\!\! \int_{\Omega_e} 6 \, \vec{u} \scall \overrightarrow{\text{grad}} \, N_i \, [\rho h]_k \, \left( N_k - \vec{\alpha} \scall \overrightarrow{\text{grad}} \, N_k \right) \mathrm{d}\Omega_e
\end{align*}

Let \( \tilde{N_k} \) be the upwind shape function, so that \( \tilde{N_k} = N_k - \vec{\alpha} \scall \overrightarrow{\text{grad}} \, N_k \) thus:

\begin{equation} \label{reynolds}
\int \!\!\! \int_{\Omega_e} \frac{\rho h^3}{\mu} \overrightarrow{\text{grad}} \, p \scall \overrightarrow{\text{grad}} \, W_i \, \mathrm{d}\Omega_e =
-\int \!\!\! \int_{\Omega_e} 6 \, \widetilde{\rho h} \; \vec{u} \scall \overrightarrow{\text{grad}} \, N_i \, \mathrm{d}\Omega_e
\end{equation}

provided that \( \widetilde{\rho h} = [\rho h]_k \, \tilde{N}_k \).

Let \( R_i \) be the residual of equation \eqref{reynolds} defined as:

\begin{equation*}
R_i = \int \!\!\! \int_{\Omega_e} \frac{\rho h^3}{\mu} \overrightarrow{\text{grad}} \, p \scall \overrightarrow{\text{grad}} \, N_i \, \mathrm{d}\Omega_e
- \int \!\!\! \int_{\Omega_e} 6 \, \widetilde{\rho h} \; \vec{u} \scall \overrightarrow{\text{grad}} \, N_i \, \mathrm{d}\Omega_e
\end{equation*}

Hence, defining \( \frac{\partial \rho}{\partial p_j} \) as \( \left[ \frac{\partial \rho}{\partial p} \right]_j N_j \) :

\begin{align}
\frac{\partial R_i}{\partial p_j} = R_{ij} =
+ & \int \!\!\! \int_{\Omega_e} \frac{\rho h^3}{\mu}                                                              \overrightarrow{\text{grad}} \, N_j \, \scall \, \overrightarrow{\text{grad}} \, N_i \, \mathrm{d}\Omega_e \nonumber \\
+ & \int \!\!\! \int_{\Omega_e} \frac{     h^3}{\mu} \left[ \frac{\partial \rho}{\partial p} \right]_j \!  N_j \, \overrightarrow{\text{grad}} \, p   \, \scall \, \overrightarrow{\text{grad}} \, N_i \, \mathrm{d}\Omega_e \nonumber \\
- & \int \!\!\! \int_{\Omega_e} 6 \, \left[ h \frac{\partial \rho}{\partial p} \right]_j \tilde{N_j} \; \vec{u} \, \scall \, \overrightarrow{\text{grad}} \, N_i \, \mathrm{d}\Omega_e
\end{align}

and \( b_i = -R_i \).


</div>
</div>
<p></p>


