**UPWINDING COEFFICIENTS**

<!-- ========== Convection-Diffusion-Reaction Equation ========== -->

<div><button class="collapsible">Convection-Diffusion-Reaction Equation</button>
<div class="content">

Zienkiewicz, "The Finite Element Method for Fluid Dynamics", 7th, p30, p50

The simplest form of Convection-Diffusion-Reaction Equation is one in which the unknown is a scalar.
It is a set of conservation laws arising from a balance of the quantity with fluxes entering and leaving a control volume:
$$
   \frac{\partial \phi}{\partial t} +\frac{\partial (U_i \phi)}{\partial x_i} -
   \frac{\partial}{\partial x_i} \left( k\frac{\partial \phi}{\partial x_i} \right) +Q =0
$$

\(U_i\) is a known velocity field and \(\phi\) is a scalar quantity being transported by this velocity.
Diffusion can also exist and \(k\) is the diffusion coefficient. The term \(Q\) represents any external
sources of the quantity \(\phi\).

</div>
</div>
<p></p>

<!-- ========== Reynolds Equation ========== -->

<div><button class="collapsible">Reynolds Equation</button>
<div class="content">

Using the fact that \(h/L\ll 1\), that the fluid is Newtonian, and that the viscosity and density are constant
through the film direction, the following Reynolds equation is obtained:
$$
   \frac{\partial }{\partial x}\left(\frac{\rho h^3}{12\mu}\frac{\partial p}{\partial x} \right) +
   \frac{\partial }{\partial y}\left(\frac{\rho h^3}{12\mu}\frac{\partial p}{\partial y} \right)= \\
   \frac{1}{2}\frac{\partial }{\partial x}\left[ \left(U_1+U_2\right)\rho h \right]+
   \frac{1}{2}\frac{\partial }{\partial y}\left[ \left(V_1+V_2\right)\rho h \right]+
              \frac{\partial }{\partial t}\left(\rho h \right)
$$

Identifying the different terms leads to:
$$
   k=\frac{\rho h^3}{\mu} \text{ ; } \phi = p \text{ ; } U_i = 6 v_i h\frac{\partial \rho}{\partial p} \text{ with } \boldsymbol{v}=
   \begin{pmatrix}
      U_2 \\
      V_2
   \end{pmatrix}
$$

<p style="text-align:center;"><img src="../media/reynolds_conf.png" alt="Reynolds conf." width="450px"/></p>

</div>
</div>
<p></p>

<!-- ========== Streamline Petrov-Galerkin weighting (SUPG) ========== -->

<div><button class="collapsible">Streamline Petrov-Galerkin weighting (SUPG)</button>
<div class="content">

In the bidimensional case, the Peclet parameter is a vector:
$$
   \boldsymbol{Pe} = \frac{1}{k} \left( \boldsymbol{U}\frac{Le}{2} \right)=
   \begin{pmatrix}
      Pe_x \\
      Pe_y
   \end{pmatrix}
$$
where \(Le\) is the element length.

The Peclet numbers of the Reynolds equation are therefore defined as:
$$
   Pe_x = \frac{6\mu(Le/2) U_2}{\rho h^2} \frac{\partial \rho}{\partial p} \\
   Pe_y = \frac{6\mu(Le/2) V_2}{\rho h^2} \frac{\partial \rho}{\partial p}
$$

The optimal upwinding coefficient is:
$$
   \alpha = \coth{Pe} -\frac{1}{Pe}
$$
with \(Pe=||\boldsymbol{Pe}||\)

<p style="text-align:center;"><img src="../media/alpha_function.png" alt="alpha function" width="350px"/></p>

</div>
</div>
<p></p>

<!-- ========== The element length \(Le\) ========== -->

<div><button class="collapsible">The element length \(Le\)</button>
<div class="content">

Using Zienkiewicz scheme:
<p style="text-align:center;"><img src="../media/Pe_h_zien.png" alt="element length z" width="350px"/></p>

which can be adapted like this:
<p style="text-align:center;"><img src="../media/Peclet_Le.png" alt="element length" width="350px"/></p>

</div>
</div>
<p></p>


