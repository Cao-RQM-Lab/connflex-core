# Theory of operation

This page is the canonical description of what connflex-core computes. Every
symbol that a user can set is listed with its parameter name in the
[table at the end](#symbols-and-parameters).

## Release approximation and beam extraction

As a self-contained polygon-to-mechanics simulator, ConnFlex automatically
extracts mechanical components given a GDS file (or a number of generic
polygons). ConnFlex targets the silicon-on-insulator (SOI) process and
expects separate polygon objects for the device, handle, and
optionally device release hole layers.

Let $D$, $H$, and $R$ be the device, handle, and release hole regions,
respectively. Isotropic expansion and erosion of region $R$ by distance $d$
are denoted by $R \oplus d$ and $R \ominus d$.

Due to finite precision in polygon processing, the processed GDS file or
polygon objects may have unintended artifacts such as hairline gaps. A healing
process applied to the device and handle layers mitigates this:

$$
L_\text{healed} = (L_\text{raw} \oplus r_h) \ominus r_h, \tag{1}
$$

where $L_\text{raw}$ is the raw layer, $L_\text{healed}$ is the healed layer,
and $r_h$ is the healing tolerance. As a performance optimization, all polygon
layers are simplified with the Douglas–Peucker algorithm with a tolerance of
$r_h$. This process contributes no additional loss of precision as features
smaller than $r_h$ were already destroyed in the prior healing step.

Release holes may be provided in a dedicated layer or within the device layer
polygon itself. In the latter case, the release holes are identified by a
single area threshold $A_\text{min}$. Holes smaller than $A_\text{min}$ are
only used in the release estimation step below and are filled. Later rigid
body calculations assume all bodies do not contain release holes.

Surviving BOX regions after release, $O$, are regions protected by both device
and handle layers, subtracted by the release hole regions, then isotropically
eroded by release distance $d$:

$$
O = \bigl((D \cap H) - R\bigr) \ominus d. \tag{2}
$$

The largest handle polygon by area is selected as the anchor for the entire
structure. ConnFlex iterates through all polygons in $D$ and $H$, discarding
all polygons without a connection to the anchor.

All features which have a smaller feature size than $w_\text{min}$ are removed
through an erosion then an expansion by $w_\text{min}/2$,

$$
B = \left(D \ominus \tfrac{w_\text{min}}{2}\right) \oplus \tfrac{w_\text{min}}{2}. \tag{3}
$$

This identifies all rigid bodies $B$. Bodies are assembled by breadth-first
search over the layer stack $B$, $O$, $H$, in which two polygons are connected
if they lie in adjacent layers and overlap. The first search originates from
the anchor and claims every polygon it reaches. This set is treated as
immovable. The search loops until all polygons are claimed.

The regions which were removed in (3) are then beam candidates $P$. Candidates
are expanded by a small amount $\epsilon$ to ensure they overlap with the
bodies they touch,

$$
P = (D - B) \oplus \epsilon. \tag{4}
$$

Beam candidates must overlap with exactly two different bodies in order to
qualify as a real beam. Beam candidates which are discarded carry neither mass
nor stiffness. Once a beam is accepted, its length $L$ is calculated from the
Euclidean distance of the overlap region centroids, its width $w$ is calculated
by dividing its geometric area by its length, and its thickness $t$ is the
device layer thickness.

!!! note "Implementation choices"
    Release holes inside the device layer are identified and filled before the
    healing of (1), so healing runs on the filled outlines rather than on every
    hole outline. Healing adds material and never removes it, so it cannot merge
    holes; it can only close a cavity whose wall has a gap narrower than
    $2 r_h$. The healed polygons are therefore scanned once more and any cavity
    at or below $A_\text{min}$ that healing closed joins the release holes.
    Holes found in the first pass keep their drawn outlines.

    The width uses the area of the unexpanded candidate $D - B$, so the
    expansion by $\epsilon$ cannot inflate it. Beam endpoints are placed at the
    overlap centroids at the device mid-plane, $z = -t_\text{device}/2$.
    Candidate parts narrower than $r_h$ are removed by an opening before
    qualification, since the body opening in (3) does not reproduce device edges
    exactly and would otherwise leave hairlines that join separate candidates.
    Candidates whose width falls outside the configured bounds are discarded
    like any other rejected candidate. Beams carry stiffness but no mass. Release
    openings, whether found as holes in the device layer or drawn on the release
    layer, are subtracted from the device-layer mass properties of a body; the filled
    polygons are kept only for connectivity and beam extraction.

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

Each layer's generic polygon is converted via constrained Delaunay
triangulation into $n$ triangles $T_n$ with edge midpoints $\mathbf{m}_{n1}$,
$\mathbf{m}_{n2}$, $\mathbf{m}_{n3}$. Then for any integrand $f$ with degree
$\leq 2$, the estimation is exact:

$$
\int f\, dA = \sum_{n} \frac{A_{T_n}}{3} \sum_{i=1}^3 f(\mathbf{m}_{ni}), \tag{11}
$$

where $A_{T_n}$ is the area of triangle $T_n$.

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

For a solid rectangular cross section with long side $a = \max(w, t)$ and
short side $b = \min(w, t)$, the Saint-Venant torsion constant is approximated
by (Young and Budynas, *Roark's Formulas for Stress and Strain*, 7th ed., 2002)

$$
J = a b^3 \left[ \frac{1}{3} - 0.21 \,\frac{b}{a} \left( 1 - \frac{b^4}{12 a^4} \right) \right]. \tag{20}
$$

This approximation approaches $a b^3 / 3$ for a thin strip and gives
$J \approx 0.141\, a^4$ for a square. ConnFlex applies an empirical
stiffening factor for short beams,

$$
J_\text{eff} = J \left( 1 + \beta \frac{a^2}{L^2} \right), \qquad \beta = 0.63 \text{ by default}. \tag{21}
$$

The coefficient $\beta$ is fitted empirically to short-beam FEA results. This
torsional correction tends to unity as $L/a$ increases. The torsional term
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

ConnFlex iterates over all $N$ identified rigid bodies (not counting the
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

## Symbols and parameters

Every symbol below is a documented input. Full definitions, units, defaults,
and validation rules are on the [parameters page](parameters.md).

| Symbol | Meaning | Parameter |
| --- | --- | --- |
| $r_h$ | healing tolerance and simplification tolerance | `ExtractionParams.heal_distance`, `ExtractionParams.simplify_tolerance` |
| $A_\text{min}$ | area threshold for release holes found inside the device layer | `ExtractionParams.release_hole_max_area` |
| $d$ | release undercut distance | `ExtractionParams.release_distance` |
| $w_\text{min}$ | minimum feature width of a rigid body; also the upper width bound of a beam | `ExtractionParams.body_min_width`, `ExtractionParams.beam_width_max` |
| $\epsilon$ | beam candidate expansion | `ExtractionParams.beam_overlap_margin` |
| $t_l$ | layer thickness | `Layer.thickness` on each `Stack` layer |
| $\rho_l$ | layer density | `Material.density` |
| $E$, $G$ | Young's and shear modulus of the device layer | `Material.youngs_modulus`, `Material.shear_modulus` or `Material.poisson_ratio` |
| $J$ | torsion constant model | `Material.torsion_model` |
| $\beta$ | short-beam torsional stiffening coefficient | `Material.torsion_beta` |
