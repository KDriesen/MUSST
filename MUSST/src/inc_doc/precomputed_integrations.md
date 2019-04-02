**PRECOMPUTED INTEGRATION PARTS**

<!-- ========== Coefficients for the transformation from (\(\xi\), \(\eta\)) to \((x,y)\) ========== -->

<div><button class="collapsible">Coefficients for the transformation from (\(\xi\), \(\eta\)) to \((x,y)\)</button>
<div class="content">

\[
   x(\xi,\eta) = \sum\limits_{nodes \; i} x_i \times N_{i}(\xi,\eta) \\
   y(\xi,\eta) = \sum\limits_{nodes \; i} y_i \times N_{i}(\xi,\eta)
\]

Let \(c_i\) be defined by:
\[
   c_1 = \frac{\partial x}{\partial \xi } \; ; \;
   c_2 = \frac{\partial y}{\partial \xi } \; ; \;
   c_3 = \frac{\partial x}{\partial \eta} \; ; \;
   c_4 = \frac{\partial y}{\partial \eta} \\
   c_5 = det \begin{pmatrix}
              c_1 & c_2 \\
              c_3 & c_4
            \end{pmatrix} = c_1 c_4 - c_2 c_3
\]

It must be remembered that the functions are first order, *ie* the \(x\) derivative with respect to \(\xi\) only depends on \(\eta\).

</div>
</div>
<p></p>

<!-- ========== Computation of the coefficients for the Reynolds equation ========== -->

<div><button class="collapsible">Computation of the coefficients for the Reynolds equation</button>
<div class="content">

\[
   h(x,y)^3 = \sum\limits_{nodes \; k} h^3_{k} \times N_{k}(x,y)
\]

The same goes for \(\rho\) and \(\mu\). Therefore, at Gauss point \((i,j)\), \(\rho h^3 /\mu\) can be precalculated:
\[
\begin{align*}
   vcal(2, i, j) &= \left( \sum\limits_{nodes \; k} h^3_{k}           \times N_{k}(x_i, y_j) \right).
                    \left( \sum\limits_{nodes \; k} \rho_{k}          \times N_{k}(x_i, y_j) \right).
                    \left( \sum\limits_{nodes \; k} \frac{1}{\mu_{k}} \times N_{k}(x_i, y_j) \right)\\
                 &= (h^3)_{ij} \; (\rho)_{ij} \; \left(\frac{1}{\mu}\right)_{ij}
\end{align*}
\]

Likewise:
\[
\begin{align*}
   vcal( 3, i, j) &= 6\, (V\!x)_{ij} \\
   vcal( 4, i, j) &= 6\, (V\!y)_{ij} \\
   vcal( 5, i, j) &= -\rho_k h_k   \,.\, \tilde{N}_{k}(x_i, y_j) \\
   vcal( 6, i, j) &= K_{ij}\left[ (\rho)_{ij}  (h^3)_{ij} \left(\frac{1}{\mu}\right)_{ij} \left(\frac{\partial p}{\partial x}\right)_{ij} \right] \\
   vcal( 7, i, j) &= K_{ij}\left[ (\rho)_{ij}  (h^3)_{ij} \left(\frac{1}{\mu}\right)_{ij} \left(\frac{\partial p}{\partial y}\right)_{ij} \right] \\
   vcal( 8, i, j) &= K_{ij}\left[ 6 \; (V\!x)_{ij} \left(\rho_k h_k \tilde{N}_k (x_i, y_j)\right) \right] \\
   vcal( 9, i, j) &= K_{ij}\left[ 6 \; (V\!y)_{ij} \left(\rho_k h_k \tilde{N}_k (x_i, y_j)\right) \right] \\
   vcal(10, i, j) &= K_{ij}\left[ (h^3)_{ij} \left(\frac{1}{\mu}\right)_{ij} \left(\frac{\partial p}{\partial x}\right)_{ij} \right] \\
   vcal(11, i, j) &= K_{ij}\left[ (h^3)_{ij} \left(\frac{1}{\mu}\right)_{ij} \left(\frac{\partial p}{\partial y}\right)_{ij} \right] \\
   vcal(12, i, j) &= K_{ij}\left[ (p)_{ij} \right]\\
   vcal(13, i, j) &= K_{ij}\left[ -\frac{1}{2} (h)_{ij} \left(\frac{\partial p}{\partial x}\right)_{ij} -(\mu)_{ij} (V\!x)_{ij} \left( \frac{1}{h} \right)_{ij} \right]
\end{align*}
\]

with:
\[
   vcal(1, i, j) = K_{ij} = c_5  w_i w_j \text{ ; } w_i \text{ Gauss ith weighting coefficient}
\]

</div>
</div>
<p></p>
