# Mass, stiffness, and modes

With the bodies and beams known, the mechanical model is a set of rigid bodies,
each with six degrees of freedom about its center of mass, joined by
Euler-Bernoulli beams welded to them. This page gives the mass and inertia of a
body from its layers, the stiffness that a beam presents to the two bodies it
joins, the assembled system, and the eigenproblem that gives its resonant
modes. `extract` computes the mass properties of every body, `Structure.modal`
assembles the system and solves for its modes, and `Structure.static` solves
it under constant loads.

## Body mass and moment of inertia

The top of the device layer is taken to be $z = 0$, with layer thicknesses
$t_l$. A body may consist of multiple layers with different geometries, with
the total mass of each layer $l$ given by

$$
m_l = \rho_l t_l A_l, \tag{5}
$$

where $\rho_l$ is the density of the layer material, and $A_l$ is the layer
polygon area. The total mass of the body $m$ is then

$$
m = \sum_l m_l, \tag{6}
$$

and the center of mass $\mathbf{c}$ of the body is

$$
\mathbf{c} = \frac{1}{m} \sum_l m_l \begin{bmatrix} \bar{x}_l \\ \bar{y}_l \\ z_l \end{bmatrix}, \tag{7}
$$

where $\bar{x}_l$ and $\bar{y}_l$ are the in-plane 2D centroid, and $z_l$ is
the mid-plane height of each layer, given by

$$
z_l = -\frac{t_l}{2} - \sum_{i < l} t_i. \tag{8}
$$

The moments of inertia in three dimensions of a body about its center of mass
are given by

$$
\begin{aligned}
I_{xx} &= \int \rho \,(y^2 + z^2)\, dV, &
I_{xy} &= -\int \rho \, x y \, dV, \\
I_{yy} &= \int \rho \,(x^2 + z^2)\, dV, &
I_{xz} &= -\int \rho \, x z \, dV, \\
I_{zz} &= \int \rho \,(x^2 + y^2)\, dV, &
I_{yz} &= -\int \rho \, y z \, dV,
\end{aligned} \tag{9}
$$

with all coordinates shifted to the body center of mass. Since every layer is
assumed to have vertical sidewalls along $z$ and is of uniform material, the
integrals can be analytically decoupled to

$$
\begin{aligned}
I_{lxx} &= \rho_l \left[ (z_{lt} - z_{lb}) \int y^2\, dA + \frac{z_{lt}^3 - z_{lb}^3}{3} \int dA \right], \\
I_{lyy} &= \rho_l \left[ (z_{lt} - z_{lb}) \int x^2\, dA + \frac{z_{lt}^3 - z_{lb}^3}{3} \int dA \right], \\
I_{lzz} &= \rho_l \,(z_{lt} - z_{lb}) \left( \int x^2\, dA + \int y^2\, dA \right), \\
I_{lxy} &= -\rho_l \,(z_{lt} - z_{lb}) \int x y\, dA, \\
I_{lxz} &= -\rho_l \,\frac{z_{lt}^2 - z_{lb}^2}{2} \int x\, dA, \\
I_{lyz} &= -\rho_l \,\frac{z_{lt}^2 - z_{lb}^2}{2} \int y\, dA,
\end{aligned} \tag{10}
$$

where $z_{lt}$ and $z_{lb}$ are the top and bottom z-heights of layer $l$ with
respect to the body center of mass.

Every integral in (5), (7), and (10) is a polynomial of degree at most two
over a region bounded by polygons, so it has a closed form in the vertices.
Let a ring have vertices $(x_i, y_i)$ for $i = 0, \dots, n - 1$, closed by
$(x_n, y_n) = (x_0, y_0)$ and traversed counterclockwise, and let

$$
a_i = x_i y_{i+1} - x_{i+1} y_i \tag{11a}
$$

be twice the signed area of the triangle spanned by the origin and edge $i$.
By Green's theorem, the integrals over the region enclosed by the ring are

$$
\begin{aligned}
\int dA &= \frac{1}{2} \sum_{i=0}^{n-1} a_i, \\
\int x\, dA &= \frac{1}{6} \sum_{i=0}^{n-1} (x_i + x_{i+1})\, a_i, \\
\int y\, dA &= \frac{1}{6} \sum_{i=0}^{n-1} (y_i + y_{i+1})\, a_i, \\
\int x^2\, dA &= \frac{1}{12} \sum_{i=0}^{n-1} \left(x_i^2 + x_i x_{i+1} + x_{i+1}^2\right) a_i, \\
\int y^2\, dA &= \frac{1}{12} \sum_{i=0}^{n-1} \left(y_i^2 + y_i y_{i+1} + y_{i+1}^2\right) a_i, \\
\int x y\, dA &= \frac{1}{24} \sum_{i=0}^{n-1} \left(2 x_i y_i + x_i y_{i+1} + x_{i+1} y_i + 2 x_{i+1} y_{i+1}\right) a_i.
\end{aligned} \tag{11b}
$$

Each term of (11b) is the exact integral of its integrand over the triangle
spanned by the origin and edge $i$, and these triangles cover the enclosed
region exactly once with sign, so the sums are exact. Traversing the ring
clockwise negates every sum; the sums of each ring are therefore multiplied
by the sign of its own $\tfrac{1}{2} \sum_i a_i$, so the drawn orientation
does not matter. The integrals over a polygon are those of its exterior ring
minus those of each hole, and the integrals over a layer region are the sum
over its polygons. $A_l$ in (5) is $\int dA$, and the centroid in (7) is
$\bar{x}_l = \int x\, dA / A_l$ and $\bar{y}_l = \int y\, dA / A_l$. The
integrals of (10) about the center of mass follow from those about the origin
by the shift to a point $\mathbf{c} = (c_x, c_y)$,

$$
\begin{aligned}
\int (x - c_x)\, dA &= \int x\, dA - c_x \int dA, \\
\int (x - c_x)^2\, dA &= \int x^2\, dA - 2 c_x \int x\, dA + c_x^2 \int dA, \\
\int (x - c_x)(y - c_y)\, dA &= \int x y\, dA - c_x \int y\, dA - c_y \int x\, dA + c_x c_y \int dA,
\end{aligned} \tag{11c}
$$

and likewise for $y$.

The three dimensional degrees of freedom of a body are defined as

$$
\mathbf{q}_b = \begin{bmatrix} u_x & u_y & u_z & \theta_x & \theta_y & \theta_z \end{bmatrix}^\top, \tag{12}
$$

where $u_x$, $u_y$, and $u_z$ are the displacement of the body's center of
mass in the global frame along each of the axes $(x, y, z)$, and $\theta_x$,
$\theta_y$, and $\theta_z$ are the rotation of the body about the axis through
the body's center of mass and parallel to each of the axes $(x, y, z)$. The
associated forces are then

$$
\mathbf{f}_b = \begin{bmatrix} f_x & f_y & f_z & \tau_x & \tau_y & \tau_z \end{bmatrix}^\top, \tag{13}
$$

where $f_x$, $f_y$, and $f_z$ are forces along each of the axes $(x, y, z)$,
and $\tau_x$, $\tau_y$, and $\tau_z$ are the torques about the axis through
the body's center of mass and parallel to each of the axes $(x, y, z)$.

With this convention, the mass and moment of inertia matrix of a body $b$ is a
$6 \times 6$ matrix:

$$
\mathbf{M}_b = \begin{bmatrix} m \mathbf{1}_3 & \mathbf{0} \\ \mathbf{0} & \mathbf{I}_b \end{bmatrix}, \tag{14}
$$

where $\mathbf{1}_3$ is the $3 \times 3$ identity matrix, and $\mathbf{I}_b$
is the body's moment of inertia matrix

$$
\mathbf{I}_b = \begin{bmatrix} I_{xx} & I_{xy} & I_{xz} \\ I_{xy} & I_{yy} & I_{yz} \\ I_{xz} & I_{yz} & I_{zz} \end{bmatrix}. \tag{15}
$$

## Beam stiffness

For a generic beam of length $L$, in-plane width $w$, and thickness $t$,
define the local beam x-axis to be along $L$, local y-axis to be in-plane, and
local z-axis to be out-of-plane. With DOF order defined in (12), the classical
Euler-Bernoulli $12 \times 12$ stiffness matrix in the beam's local frame is

$$
\mathbf{K}^\text{loc} = \begin{bmatrix} \mathbf{k}_{00} & \mathbf{k}_{01} \\ \mathbf{k}_{01}^\top & \mathbf{k}_{11} \end{bmatrix}, \tag{16}
$$

$$
\mathbf{k}_{00} = \begin{bmatrix}
k_a & 0 & 0 & 0 & 0 & 0 \\
0 & k_y & 0 & 0 & 0 & c_z \\
0 & 0 & k_z & 0 & -c_y & 0 \\
0 & 0 & 0 & k_t & 0 & 0 \\
0 & 0 & -c_y & 0 & g_y & 0 \\
0 & c_z & 0 & 0 & 0 & g_z
\end{bmatrix},
\qquad
\mathbf{k}_{11} = \begin{bmatrix}
k_a & 0 & 0 & 0 & 0 & 0 \\
0 & k_y & 0 & 0 & 0 & -c_z \\
0 & 0 & k_z & 0 & c_y & 0 \\
0 & 0 & 0 & k_t & 0 & 0 \\
0 & 0 & c_y & 0 & g_y & 0 \\
0 & -c_z & 0 & 0 & 0 & g_z
\end{bmatrix}, \tag{17}
$$

$$
\mathbf{k}_{01} = \begin{bmatrix}
-k_a & 0 & 0 & 0 & 0 & 0 \\
0 & -k_y & 0 & 0 & 0 & c_z \\
0 & 0 & -k_z & 0 & -c_y & 0 \\
0 & 0 & 0 & -k_t & 0 & 0 \\
0 & 0 & c_y & 0 & g_y/2 & 0 \\
0 & -c_z & 0 & 0 & 0 & g_z/2
\end{bmatrix}, \tag{18}
$$

$$
\begin{aligned}
k_a &= \frac{E w t}{L}, & k_t &= \frac{G J}{L}, \\
k_y &= \frac{E t w^3}{L^3}, & c_z &= \frac{E t w^3}{2 L^2}, & g_z &= \frac{E t w^3}{3 L}, \\
k_z &= \frac{E w t^3}{L^3}, & c_y &= \frac{E w t^3}{2 L^2}, & g_y &= \frac{E w t^3}{3 L},
\end{aligned} \tag{19}
$$

where $E$ is the Young's modulus of the beam material, $G$ is the shear
modulus, and $J$ is the torsion constant.

The element is that of Euler-Bernoulli theory: plane sections stay plane and
square to the axis, so shear deformation is left out. Against a Timoshenko
beam every bending term of (19) is too large by $1 + \Phi$ with $\Phi = 12 E I
/ (\kappa G A L^2)$: out of plane $\Phi \approx 3.07\, (t / L)^2$ for silicon
with $\nu = 0.28$ and $\kappa = 5/6$, so a beam eight thicknesses long is 5
percent too stiff out of plane and one four thicknesses long 19 percent, and
in the plane the same with $w$ for $t$. Bending uses $E$ as given in both
directions: in-plane bending of a beam much deeper than wide approaches the
plate modulus $E / (1 - \nu^2)$, 8.5 percent higher at $\nu = 0.28$, and one
isotropic $E$ serves beams of every orientation, where silicon's differs by
up to 23 percent between $\langle 110 \rangle$ and $\langle 100 \rangle$.
Beams carry no mass, no residual or thermal stress, and no root or anchor
compliance: their ends are welded to rigid bodies at the body face, and the
anchor does not move.

For a solid rectangular cross section with long side $a = \max(w, t)$ and
short side $b = \min(w, t)$, the Saint-Venant torsion constant is approximated
by (Young and Budynas, *Roark's Formulas for Stress and Strain*, 7th ed., 2002)

$$
J = a b^3 \left[ \frac{1}{3} - 0.21 \,\frac{b}{a} \left( 1 - \frac{b^4}{12 a^4} \right) \right]. \tag{20}
$$

This approximation approaches $a b^3 / 3$ for a thin strip and gives
$J \approx 0.141\, a^4$ for a square. connflex-core applies an empirical
stiffening factor for short beams,

$$
J_\text{eff} = J \left( 1 + \beta \frac{a^2}{L^2} \right), \qquad \beta = 0.63 \text{ by default}. \tag{21}
$$

The coefficient $\beta$ has no source here: it was set from short-beam
finite-element results whose data and range of $L / a$ were not kept, no
test checks it beyond its own algebra, and restrained warping, the mechanism
behind such stiffening, scales as $a / L$ for a thin strip rather than
$a^2 / L^2$. Take the correction as unverified. It tends to unity as $L / a$
increases. The torsional term
$k_t$ in (19) uses $J_\text{eff}$, computed per beam from its own $w$, $t$,
and $L$. Setting $\beta = 0$ disables the correction.

Take a generic beam connected to bodies $A$ and $B$ at points $\mathbf{p}_0$
and $\mathbf{p}_1$ (global frame), respectively. Since it is assumed that a
beam only ever lives in the device layer, its rotation relative to the global
frame is always only a single $\alpha$. Hence the rotation transformation
from the beam's local frame to the body (global) frame is

$$
\mathbf{R}_3 = \begin{bmatrix} \cos \alpha & -\sin \alpha & 0 \\ \sin \alpha & \cos \alpha & 0 \\ 0 & 0 & 1 \end{bmatrix}, \tag{22}
$$

where

$$
\alpha = \operatorname{atan2}(\mathbf{v} \cdot \hat{\jmath}, \mathbf{v} \cdot \hat{\imath}), \qquad \mathbf{v} = \mathbf{p}_1 - \mathbf{p}_0. \tag{23}
$$

$\mathbf{R}_3$ is orthogonal, therefore the same transformation is
identically applied to displacement ($\mathbf{u}$), rotation
($\boldsymbol{\theta}$), force ($\mathbf{f}$), and torque
($\boldsymbol{\tau}$), hence the rotation operation for the full DOF is

$$
\mathbf{R}_6 = \begin{bmatrix} \mathbf{R}_3 & \mathbf{0} \\ \mathbf{0} & \mathbf{R}_3 \end{bmatrix}. \tag{24}
$$

When rigid body $A$ translates by $\mathbf{u}_A$ and rotates by
$\boldsymbol{\theta}_A$, a point on this body defined by vector $\mathbf{d}_A$
away from the body's center of mass $\mathbf{c}_A$ moves by

$$
\mathbf{u}(\mathbf{c}_A + \mathbf{d}_A) = \mathbf{u}_A + \boldsymbol{\theta}_A \times \mathbf{d}_A. \tag{25}
$$

For the aforementioned beam, the attachment point vector is then

$$
\mathbf{d}_A = \mathbf{p}_0 - \mathbf{c}_A. \tag{26}
$$

Writing this relationship in terms of the body and beam DOFs, and converting
the cross product $\boldsymbol{\theta}_A \times \mathbf{d}_A$ into a
skew-symmetric matrix gives

$$
\begin{bmatrix} \mathbf{u}(\mathbf{p}_0) \\ \boldsymbol{\theta}_A \end{bmatrix}
= \begin{bmatrix} \mathbf{1}_3 & -[\mathbf{d}_A]_\times \\ \mathbf{0} & \mathbf{1}_3 \end{bmatrix}
\begin{bmatrix} \mathbf{u}_A \\ \boldsymbol{\theta}_A \end{bmatrix}
= \mathbf{S}_A^\top \mathbf{q}_A, \tag{27}
$$

where

$$
[\mathbf{d}_A]_\times = \begin{bmatrix} 0 & -d_z & d_y \\ d_z & 0 & -d_x \\ -d_y & d_x & 0 \end{bmatrix}. \tag{28}
$$

In order to transform the body displacement/rotation DOFs onto one beam end's
frame, apply a translation then a rotation:

$$
\mathbf{q}_0 = \mathbf{R}_6^\top \mathbf{S}_A^\top \mathbf{q}_A = \mathbf{T}_A^\top \mathbf{q}_A. \tag{29}
$$

By the principle of virtual work, the resultant beam end force/torque must
transform by the transpose back to the global frame:

$$
\mathbf{f}_A = \mathbf{S}_A \mathbf{R}_6 \mathbf{f}_0 = \mathbf{T}_A \mathbf{f}_0. \tag{30}
$$

The construction is identical for end 1 of the beam ($\mathbf{p}_1$,
$\mathbf{q}_1$, $\mathbf{f}_1$, $\mathbf{S}_B$, $\mathbf{T}_B$). The full
transformation matrix for a beam with two ends is then

$$
\mathbf{T} = \begin{bmatrix} \mathbf{T}_A & \mathbf{0} \\ \mathbf{0} & \mathbf{T}_B \end{bmatrix}
= \begin{bmatrix} \mathbf{S}_A \mathbf{R}_6 & \mathbf{0} \\ \mathbf{0} & \mathbf{S}_B \mathbf{R}_6 \end{bmatrix}. \tag{31}
$$

When an end of a beam is welded to a rigid body, the relationship is
represented by transforming the body's global DOFs onto the beam end's local
DOFs, applying the beam's local stiffness $\mathbf{K}^\text{loc}$, then
transforming the resulting forces back onto the body's global frame:

$$
\begin{bmatrix} \mathbf{f}_A \\ \mathbf{f}_B \end{bmatrix}
= \mathbf{T} \begin{bmatrix} \mathbf{f}_0 \\ \mathbf{f}_1 \end{bmatrix}
= \mathbf{T} \mathbf{K}^\text{loc} \begin{bmatrix} \mathbf{q}_0 \\ \mathbf{q}_1 \end{bmatrix}
= \mathbf{T} \mathbf{K}^\text{loc} \mathbf{T}^\top \begin{bmatrix} \mathbf{q}_A \\ \mathbf{q}_B \end{bmatrix}, \tag{32}
$$

therefore

$$
\mathbf{K}^\text{eff} = \mathbf{T} \mathbf{K}^\text{loc} \mathbf{T}^\top, \tag{33}
$$

where $\mathbf{K}^\text{eff}$ is the effective stiffness felt by the body due
to the beam.

## Assemble and solve

connflex-core iterates over all $N$ identified rigid bodies (not counting the
anchor body with 0 DOFs) and assembles a $6N \times 6N$ mass matrix:

$$
\mathbf{M}^\text{full} = \operatorname{blkdiag}(\mathbf{M}_1, \dots, \mathbf{M}_N). \tag{34}
$$

Similarly, a $6N \times 6N$ stiffness matrix consisting of $N^2$ total
$6 \times 6$ blocks is built:

$$
\mathbf{K}^\text{full} = \begin{bmatrix}
\mathbf{K}_{11} & \mathbf{K}_{12} & \cdots & \mathbf{K}_{1N} \\
\mathbf{K}_{12}^\top & \mathbf{K}_{22} & \cdots & \mathbf{K}_{2N} \\
\vdots & \vdots & \ddots & \vdots \\
\mathbf{K}_{1N}^\top & \mathbf{K}_{2N}^\top & \cdots & \mathbf{K}_{NN}
\end{bmatrix}, \tag{35}
$$

where each block $\mathbf{K}_{ij}$ is the sum of all contributing blocks of
$\mathbf{K}^\text{eff}$ from every beam which connects to body $i$ and $j$.
The anchor body has all associated blocks discarded. For example, a beam
connected to body $a$ and body $b$ will additively contribute 4 blocks:

$$
\begin{bmatrix} \mathbf{K}_{aa} & \mathbf{K}_{ab} \\ \mathbf{K}_{ab}^\top & \mathbf{K}_{bb} \end{bmatrix}
= \begin{bmatrix}
\mathbf{T}_a \mathbf{k}_{00}^\text{loc} \mathbf{T}_a^\top & \mathbf{T}_a \mathbf{k}_{01}^\text{loc} \mathbf{T}_b^\top \\
\mathbf{T}_b \mathbf{k}_{01}^{\text{loc}\top} \mathbf{T}_a^\top & \mathbf{T}_b \mathbf{k}_{11}^\text{loc} \mathbf{T}_b^\top
\end{bmatrix}. \tag{36}
$$

If body $a$ is the anchor body, the contribution is then

$$
\begin{bmatrix} \mathbf{0} & \mathbf{0} \\ \mathbf{0} & \mathbf{K}_{bb} \end{bmatrix}. \tag{37}
$$

Finally, the resonances of this reduced-order system are found by solving the
eigenvalue problem

$$
\mathbf{K}^\text{full} \boldsymbol{\phi}_n = \omega_n^2 \mathbf{M}^\text{full} \boldsymbol{\phi}_n, \qquad n = 1, \dots, 6N, \tag{38}
$$

where $\omega_n$ is the resonant angular frequency and $\boldsymbol{\phi}_n$
is the mechanical mode.

!!! warning "Beams stay linear"
    $\mathbf{K}^\text{full}$ is the stiffness of small deflection, and every
    analysis takes it at any displacement: static equilibrium, small signal
    about it, and transients. The field follows the displaced bodies
    exactly; the beams do not. A beam held at both ends stretches as it
    bends, and stiffens once its deflection is no longer small against its
    width in the plane of bending; any beam stiffens or softens once its
    deflection or the turn of its ends is no longer small against its length
    or a radian, and the direction it pulls in turns with the bodies. None
    of that is in the model, and nothing warns of it. Compare the travel of
    a result with the range of its flexures before trusting it, most of all
    for motors driven over many pitches. Shear deformation, beam mass,
    residual stress, thermal expansion, root and anchor compliance, and the
    anisotropy of silicon are left out as well; the beam stiffness section
    says what each costs.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $t_l$, $\rho_l$ | thickness and density of layer $l$ | `Layer.thickness`, `Material.density` |
| $m$, $\mathbf{c}$, $\mathbf{I}_b$ | mass, center of mass, and inertia tensor of a body | `Body.inertial.mass`, `Body.inertial.center`, `Body.inertial.inertia` |
| $\mathbf{M}_b$ | 6 by 6 mass matrix of a body | `Body.inertial.matrix` |
| $\mathbf{q}_b$, $\mathbf{f}_b$ | degrees of freedom and generalized forces of a body | rows of `ModalResult.dofs`; loads of `Structure.static` |
| $E$, $G$, $\nu$ | Young's modulus, shear modulus, Poisson's ratio of the device layer | `Material.youngs_modulus`, `Material.shear_modulus`, `Material.poisson_ratio` |
| $J$, $J_\text{eff}$, $\beta$ | torsion constant, its short-beam correction, and the correction's coefficient | `Material.torsion_beta` |
| $\mathbf{K}^\text{loc}$, $\mathbf{K}^\text{eff}$ | stiffness of a beam in its own frame and as the bodies feel it | `mechanics.stiffness.local_stiffness`, `effective_stiffness` |
| $\mathbf{M}^\text{full}$, $\mathbf{K}^\text{full}$ | assembled mass and stiffness, anchor removed | `ModalResult.mass`, `ModalResult.stiffness` |
| $\omega_n$, $\boldsymbol{\phi}_n$ | resonant angular frequency and mode shape | `ModalResult.frequencies_hz` as $\omega_n / 2\pi$, columns of `ModalResult.modes` |
