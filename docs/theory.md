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
qualify as a real beam. A discarded candidate that overlaps exactly one body is
a stub, such as a comb finger: it is rigidly attached to that body and
contributes mass and moment of inertia, but no stiffness. Discarded
candidates that overlap no body or three or more bodies carry neither mass
nor stiffness. Once a beam is accepted, its length $L$ is calculated from the
Euclidean distance of the overlap region centroids, its width $w$ is calculated
by dividing its geometric area by its length, and its thickness $t$ is the
device layer thickness.

!!! note "Implementation choices"
    Every polygon is put in canonical vertex order when a layout is built: the
    exterior ring runs clockwise and starts at its lexicographically smallest
    vertex, each hole runs counterclockwise and starts at its own smallest
    vertex, and holes are listed in a fixed order given by their coordinates.
    Polygons keep the order of the input. The Douglas–Peucker simplification
    that follows (1) treats a ring as an open sequence from its start vertex
    around to the same vertex: it keeps both ends, keeps the vertex farthest
    from the chord between them when that distance exceeds $r_h$, and recurses
    on the two halves. The retained vertices therefore depend on the start
    vertex, and the same outline listed from another vertex simplifies to a
    different polygon within $r_h$ of the first. Erosion and expansion do not
    depend on where a ring starts or which way it runs, so with canonical
    order every layer downstream of (1) is a function of the drawn shapes
    alone.

    Release holes inside the device layer are identified and filled before the
    healing of (1), so healing runs on the filled outlines rather than on every
    hole outline. Healing adds material and never removes it, so it cannot merge
    holes; it can only close a cavity whose wall has a gap narrower than
    $2 r_h$. The healed polygons are therefore scanned once more and any cavity
    at or below $A_\text{min}$ that healing closed joins the release holes.
    Holes found in the first pass keep their drawn outlines. Interior rings of a
    valid polygon are disjoint, and a cavity closed by healing lay in open
    space before it, so the release holes never overlap one another, and a
    drawn release layer is unioned by its reader; a union runs only when both
    are present. The resulting openings are built
    once and serve both the oxide estimate in (2) and the mass properties.
    Equation (2) is evaluated as $((D \cap H) \ominus d) - (R \oplus d)$, which is
    the same set: a point survives the erosion of a difference exactly when its
    disk lies in the minuend and misses the subtrahend. The protected region is
    eroded on its own, and only the openings within $d$ of the eroded region are
    dilated and subtracted, rather than eroding one polygon that carries every
    opening as a hole.

    The width uses the area of the unexpanded candidate $D - B$, so the
    expansion by $\epsilon$ cannot inflate it. Beam endpoints are placed at the
    overlap centroids at the device mid-plane, $z = -t_\text{device}/2$.
    Candidate parts narrower than $r_h$ are removed by an opening before
    qualification, since the body opening in (3) does not reproduce device edges
    exactly and would otherwise leave hairlines that join separate candidates.
    Candidates whose width falls outside the configured bounds are discarded
    like any other rejected candidate. Beams carry stiffness but no mass.

    The opening in (3) and the sliver removal on the candidates are applied
    polygon by polygon. An opening never leaves its polygon, so the parts of a
    disjoint set stay disjoint and the result is the same as opening the set as
    a whole, without the union that a single buffer over the set performs.

    Stubs are kept apart from the body polygons. A body's device polygons are
    the output of the opening in (3): they define connectivity through the
    stack, the candidate regions, and the net tags, and widening them with
    stubs would change what qualifies as a beam. A stub is therefore stored on
    its body as a separate polygon set and still appears among the rejected
    candidates with its reason, since it is not a beam. For mass, center of
    mass, and inertia the body is integrated over its device polygons and its
    stubs alike, with release openings subtracted from both. Candidates that
    touch three or more bodies have no single owner and candidates between two
    bodies that fail the width bounds are neither body nor beam, so both stay
    massless. No parameter changes this: a massless comb is wrong, not a
    modeling option. Release
    openings, whether found as holes in the device layer or drawn on the release
    layer, are subtracted from the device-layer mass properties of a body; the filled
    polygons are kept only for connectivity and beam extraction. The subtraction
    works on the moments: the sums of (11b) are evaluated for the filled
    device polygons of a body and for the openings it contains, and the latter
    are taken away. An opening that crosses the body boundary is clipped to it
    first. The result equals integrating over the perforated polygon, and all
    rings of a geometry are summed in one vectorized pass, so a layout with a
    hundred thousand openings costs one pass over their vertices and no
    polygon is triangulated.

    The sums of (11b) are evaluated about the center $\mathbf{c}_0$ of the
    bounding box of the region and moved to the origin by (11c) with the shift
    reversed. About the global origin every term of (11b) is of order
    $|\mathbf{c}_0|^2 A$, whereas the central moments of (10) are of order
    $\sigma^2 A$ for a body of extent $\sigma$, so a body a few millimeters
    from the origin would lose about seven significant digits of inertia in
    the sum. Summed about $\mathbf{c}_0$, the only cancellation left is the
    shift to the center of mass in (11c).

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

## Electrodes and links

The electrical model is a network whose nodes are electrodes and whose
branches are links. Its conductors are the connected regions of the healed
device and handle layers that survive pruning, and nets are the groups of
conductors joined through vias.

Let $B_j$ be the device polygons of the bodies from (3) and $P_k$ the beam
candidates $D - B$ of (4). Candidate $P_k$ touches body polygon $B_j$ when
its expansion overlaps it,

$$
P_k \sim B_j \iff \operatorname{area}\bigl((P_k \oplus \epsilon) \cap B_j\bigr) > 0, \tag{39}
$$

which is the test that qualifies beams. Let $T_k = \{\, j : P_k \sim B_j \,\}$
be the set of body polygons that candidate $k$ touches. The electrode of body
polygon $B_j$ is that polygon together with every candidate that touches it
and no other body polygon,

$$
E_j = B_j \cup \bigcup_{k \,:\, T_k = \{j\}} P_k. \tag{40}
$$

Every handle polygon is an electrode by itself. A candidate with
$T_k = \emptyset$, such as a narrow trace that stands on its own oxide, is an
electrode by itself as well. An electrode moves with one body and lies in one
conductor. For $E_j$ these are the body that owns $B_j$ and the conductor that
contains it. A handle electrode moves with the body that owns its polygon. A
candidate that touches no body polygon moves with the body whose oxide it
overlaps over the largest area.

Links conduct between electrodes:

- A beam joins two electrodes, one on each of its two bodies. Where it touches
  several polygons of one body, it joins the largest.
- A rejected candidate with $|T_k| \geq 2$ is a bridge. It joins the electrode
  of the largest polygon it touches to the electrode of each other polygon in
  $T_k$.
- A via joins device electrodes to handle electrodes. It reaches one electrode
  of each device conductor it overlaps: the one it overlaps over the largest
  area, or the nearest one when the via lies wholly on a beam or a bridge. It
  joins each of these to the electrode of every handle polygon it overlaps,
  with one link per pair. A via that overlaps no conductor on one side joins
  nothing.

A beam conducts along its length through its full cross-section, so its link
has resistance

$$
R = \frac{\rho_e L}{w t}, \tag{41}
$$

where $\rho_e$ is the electrical resistivity of the device layer and $L$, $w$,
and $t$ are the length, width, and thickness of the beam as extracted: $L$ is
the distance between the overlap centroids and $w$ is the area $A$ of the
unexpanded candidate divided by $L$, so $R = \rho_e L^2 / (A t)$. Bridges are
ideal, $R = 0$. Each link of a via has resistance $R_v$, one value supplied
for all vias and independent of via area, which stands for both contacts and
the plug between them. It is zero unless supplied, which makes vias ideal.

Resistance inside an electrode is neglected, so every electrode is an
equipotential. A beam contributes $L / w$ squares of the sheet resistance
$\rho_e / t$, typically tens to hundreds, whereas current spreading from a
beam root across a plate contributes of the order of one square. The same
neglect covers stubs: a comb finger of length $\ell$ and width $w_f$ adds up
to $\ell / w_f$ squares between its tip and its root, which is not modeled.

Links never join electrodes of different nets: a beam or a bridge lies within
one conductor, and a via joins conductors of one net. Conversely, the
electrodes of one net are joined through links, because the pieces of a
conductor that the opening in (3) separates are joined by the candidates
between them. There are two exceptions. A candidate narrower than $r_h$ is
removed as a sliver before qualification and leaves no link, and a beam whose
end touches several polygons of one body joins only the largest.

!!! note "Implementation choices"
    Neither a body nor a net can serve as a node. The anchor body carries every
    bonded pad and every stator, which sit on different nets, and one net runs
    through every body that beams join to its pad. An electrode is the largest
    piece of silicon that has both one potential and one rigid motion.

    A candidate that touches two polygons of the same body is a stub of that
    body for mass, since it moves with it, and a bridge for conduction, since
    it joins two electrodes. The silicon of a bridge or a beam belongs to no
    electrode.

    Device electrodes are listed first and handle electrodes after them, each
    ordered by body and then by descending area of the body polygon.

    A via carries one resistance rather than a value per area because the
    `vias` layer states only that device and handle are connected, not how the
    contact is made. A via that straps several device electrodes joins each of
    them to the handle through $R_v$, so two of them see $2 R_v$ between each
    other.

    Resistivity of an SOI device layer spans orders of magnitude with doping,
    so it has no default. Without it the resistance of a beam link is left
    unknown rather than guessed, and extraction is otherwise unchanged.

## Plates between device and handle

Device silicon that lies over handle silicon forms a capacitor with it:
through buried oxide where oxide survives, and across the gap that release
leaves elsewhere. Silicon is treated as a perfect conductor, so an electrode
has one potential, and a plate is a pair of a device electrode $E_a$ and a
handle electrode $E_b$ whose regions overlap. With $O$ the surviving oxide of
(2) and $R$ the release openings, the overlap and its bonded and released
parts are

$$
U_{ab} = E_a \cap E_b, \qquad O_{ab} = U_{ab} \cap O, \qquad G_{ab} = U_{ab} - O - R. \tag{42}
$$

Oxide and gap both have the thickness of the box layer, $g_0 = t_\text{box}$.
The field between the two faces is taken as uniform and vertical, so at rest

$$
C_{ab} = C_\text{ox} + \frac{\varepsilon A}{g_0}, \qquad
C_\text{ox} = \frac{\varepsilon_0\, \varepsilon_\text{ox} \operatorname{area}(O_{ab})}{g_0}, \qquad
\varepsilon = \varepsilon_0\, \varepsilon_a, \tag{43}
$$

where $A = \operatorname{area}(G_{ab})$, $\varepsilon_\text{ox}$ is the
relative permittivity of the box material, $\varepsilon_a$ is that of the
medium in the gap, and $\varepsilon_0 = 8.8541878188 \times 10^{-12}$ F/m.
Release openings carry no charge: device silicon is absent there.

### Local coordinates

Oxide bonds the silicon on both of its faces, so $C_\text{ox}$ does not
change with motion. The released part does. Let the device electrode move
with body $A$ and the handle electrode with body $B$, with centers of mass
$\mathbf{c}_A$ and $\mathbf{c}_B$ and degrees of freedom $\mathbf{q}_A$ and
$\mathbf{q}_B$ of (12). A small rigid motion of a body with center
$\mathbf{c}$ displaces its material at $(x, y, z)$ by
$\mathbf{u} + \boldsymbol{\theta} \times (\mathbf{r} - \mathbf{c})$,

$$
\begin{aligned}
\delta x &= u_x + \theta_y (z - c_z) - \theta_z (y - c_y), \\
\delta y &= u_y + \theta_z (x - c_x) - \theta_x (z - c_z), \\
\delta z &= u_z + \theta_x (y - c_y) - \theta_y (x - c_x).
\end{aligned} \tag{44}
$$

The gap between the faces, $g = g_0 + \delta z_A - \delta z_B$, therefore
stays a plane. About a reference point $\mathbf{r}_0 = (x_0, y_0, z_0)$, with
$(x_0, y_0)$ the center of the bounding box of $U_{ab}$ and
$z_0 = -(t_\text{device} + t_\text{box} / 2)$ the mid-height of the gap, it
is

$$
g(x, y) = g_0 + h_0 + h_x (x - x_0) + h_y (y - y_0), \tag{45}
$$

with rise $h_0$ and slopes $h_x$ and $h_y$. In the plane, device silicon
moves against handle silicon by a travel $(a_x, a_y)$ at the reference point
and a rotation $a_\psi$ about the vertical through it. These six local
coordinates are the motion of body $A$ relative to body $B$ at
$\mathbf{r}_0$, linear in the degrees of freedom of both,

$$
\mathbf{s} = \begin{bmatrix} a_x \\ a_y \\ a_\psi \\ h_0 \\ h_x \\ h_y \end{bmatrix}
= \mathbf{J}(\mathbf{c}_A)\, \mathbf{q}_A - \mathbf{J}(\mathbf{c}_B)\, \mathbf{q}_B, \qquad
\mathbf{J}(\mathbf{c}) = \begin{bmatrix}
1 & 0 & 0 & 0 & d_z & -d_y \\
0 & 1 & 0 & -d_z & 0 & d_x \\
0 & 0 & 0 & 0 & 0 & 1 \\
0 & 0 & 1 & d_y & -d_x & 0 \\
0 & 0 & 0 & 0 & -1 & 0 \\
0 & 0 & 0 & 1 & 0 & 0
\end{bmatrix}, \qquad \mathbf{d} = \mathbf{r}_0 - \mathbf{c}, \tag{46}
$$

which is (44) evaluated at $\mathbf{r}_0$ for each body. Both bodies are
evaluated at the same point, so a rigid motion common to the two gives
$\mathbf{s} = 0$.

### Gap and tilt

While the faces stay nearly parallel the field remains locally that of a
parallel plate with the local gap. With $\bar{g} = g_0 + h_0$ the gap at the
reference point, $\mathbf{h} = (h_x, h_y)^\top$, and
$\boldsymbol{\rho} = (x - x_0,\; y - y_0)^\top$, capacitance of the released
part is

$$
\varepsilon \int_{G_{ab}} \frac{dA}{g(x, y)}
= \varepsilon \int_{G_{ab}} \left[ \frac{1}{\bar{g}}
- \frac{\mathbf{h}^\top \boldsymbol{\rho}}{\bar{g}^2}
+ \frac{(\mathbf{h}^\top \boldsymbol{\rho})^2}{\bar{g}^3} - \cdots \right] dA, \tag{47}
$$

a series that converges while the faces do not touch,
$|\mathbf{h}^\top \boldsymbol{\rho}| < \bar{g}$ over $G_{ab}$. Its terms
need only the area moments of the released region about $\mathbf{r}_0$,

$$
A = \int_{G_{ab}} dA, \qquad
\mathbf{S} = \int_{G_{ab}} \boldsymbol{\rho}\, dA, \qquad
\mathbf{I} = \int_{G_{ab}} \boldsymbol{\rho} \boldsymbol{\rho}^\top dA, \tag{48}
$$

which are the integrals of (11b) shifted to $\mathbf{r}_0$ by (11c). The
model keeps (47) through second order in the slopes.

### In-plane travel

Travel in the plane changes capacitance by changing the overlap. Relative to
handle silicon, device silicon at $\boldsymbol{\rho}$ moves by
$\sum_k a_k \mathbf{w}_k$ over $k \in \{x, y, \psi\}$, with

$$
\mathbf{w}_x = \begin{bmatrix} 1 \\ 0 \end{bmatrix}, \qquad
\mathbf{w}_y = \begin{bmatrix} 0 \\ 1 \end{bmatrix}, \qquad
\mathbf{w}_\psi = \begin{bmatrix} -(y - y_0) \\ x - x_0 \end{bmatrix}. \tag{49}
$$

Seen from the device electrode, its own outline, its openings, and the gap
plane (45) stand still while the outline of the handle electrode moves by
the opposite amount. Only the pieces of that outline which lie under device
silicon bound the released region,

$$
\Gamma_{ab} = \partial E_b \cap (E_a - R), \tag{50}
$$

and with $\mathbf{n}$ the unit normal that points out of $E_b$, a piece
sweeps area at the rate $-\mathbf{n}^\top \mathbf{w}_k$ per unit length. The
integral over the released region of any function $f$ that stands still with
the device electrode therefore changes as

$$
\frac{\partial}{\partial a_k} \int_{G_{ab}} f\, dA = -\int_{\Gamma_{ab}} f\; \mathbf{n}^\top \mathbf{w}_k\, ds. \tag{51}
$$

Taking $f = 1$ and $f = \boldsymbol{\rho}$ defines the edge integrals

$$
\ell_k = \int_{\Gamma_{ab}} \mathbf{n}^\top \mathbf{w}_k\, ds, \qquad
\boldsymbol{\lambda}_k = \int_{\Gamma_{ab}} \boldsymbol{\rho}\; \mathbf{n}^\top \mathbf{w}_k\, ds, \tag{52}
$$

and gives the first two moments of (48) to first order in travel,

$$
A(\mathbf{a}) = A - \sum_k a_k\, \ell_k, \qquad
\mathbf{S}(\mathbf{a}) = \mathbf{S} - \sum_k a_k\, \boldsymbol{\lambda}_k. \tag{53}
$$

For a handle edge of length $w$ that passes straight under a plate with its
normal along $x$, $\ell_x = w$: travel of the plate by $a_x$ along that
normal carries a strip $w\, a_x$ off the handle. Along a straight piece of
$\Gamma_{ab}$ from point 1 to point 2, of length $L$, both
$c = \mathbf{n}^\top \mathbf{w}_k$ and a weight
$\varphi \in \{1,\, x - x_0,\, y - y_0\}$ are linear, so every integral of
(52) is a sum over pieces of

$$
\int_1^2 c\, \varphi\, ds = \frac{L}{6} \left( 2 c_1 \varphi_1 + c_1 \varphi_2 + c_2 \varphi_1 + 2 c_2 \varphi_2 \right). \tag{54}
$$

### Capacitance and its derivatives

With (53) in the first two terms of (47), capacitance of a plate is

$$
C_{ab}(\mathbf{s}) = C_\text{ox} + \varepsilon \left[ \frac{A(\mathbf{a})}{\bar{g}}
- \frac{\mathbf{h}^\top \mathbf{S}(\mathbf{a})}{\bar{g}^2}
+ \frac{\mathbf{h}^\top \mathbf{I}\, \mathbf{h}}{\bar{g}^3} \right],
\qquad \bar{g} = g_0 + h_0. \tag{55}
$$

Its accuracy differs by direction:

- Every coefficient carries its exact dependence on the rise, so vertical
  travel without tilt is exact at any amplitude,
  $C_{ab} = C_\text{ox} + \varepsilon A / (g_0 + h_0)$.
- Tilt is carried to second order. The first term left out is
  $-\varepsilon \int (\mathbf{h}^\top \boldsymbol{\rho})^3 dA / \bar{g}^4$,
  of relative size $(|\mathbf{h}|\, \ell / \bar{g})^3$ for a plate of extent
  $\ell$; it vanishes for a region symmetric about $\mathbf{r}_0$, which
  leaves an error of fourth order.
- In-plane travel is carried to first order, together with its products
  with rise and slopes: by (51) with $f = 1 / g$,
  $\partial C_{ab} / \partial a_k = -\varepsilon \int_{\Gamma_{ab}}
  \mathbf{n}^\top \mathbf{w}_k\, (1 / \bar{g} - \mathbf{h}^\top
  \boldsymbol{\rho} / \bar{g}^2)\, ds$, the handle outline under the plate
  weighted by the inverse of the local gap. Left out are curvature of the
  overlap in travel, which depends on the angles at which the two outlines
  cross and on rotation beyond first order, where (44) is itself linearized;
  the change of $\mathbf{I}$ with travel, of third order; and the product of
  travel with tilt of body $B$, since the gap plane is taken to ride with
  the device electrode. The last vanishes when body $B$ is the anchor.
- Equation (53) holds until a corner of either outline or an opening
  reaches $\Gamma_{ab}$. A release hole drawn next to a handle edge ends
  that range at once in one direction. Where a piece of handle outline runs
  along the outline of device silicon or of an opening, overlap has a kink
  and no derivative; such a piece is counted in $\Gamma_{ab}$ along the
  outline of the device electrode and not along an opening.

Gradient and Hessian are the derivatives of (55) in closed form,

$$
\begin{aligned}
\frac{\partial C_{ab}}{\partial a_k} &= \varepsilon \left[ -\frac{\ell_k}{\bar{g}}
+ \frac{\mathbf{h}^\top \boldsymbol{\lambda}_k}{\bar{g}^2} \right], \\
\frac{\partial C_{ab}}{\partial h_0} &= \varepsilon \left[ -\frac{A(\mathbf{a})}{\bar{g}^2}
+ \frac{2\, \mathbf{h}^\top \mathbf{S}(\mathbf{a})}{\bar{g}^3}
- \frac{3\, \mathbf{h}^\top \mathbf{I}\, \mathbf{h}}{\bar{g}^4} \right], \\
\frac{\partial C_{ab}}{\partial \mathbf{h}} &= \varepsilon \left[ -\frac{\mathbf{S}(\mathbf{a})}{\bar{g}^2}
+ \frac{2\, \mathbf{I}\, \mathbf{h}}{\bar{g}^3} \right],
\end{aligned} \tag{56}
$$

$$
\begin{aligned}
\frac{\partial^2 C_{ab}}{\partial a_k\, \partial a_l} &= 0, &
\frac{\partial^2 C_{ab}}{\partial h_0^2} &= \varepsilon \left[ \frac{2 A(\mathbf{a})}{\bar{g}^3}
- \frac{6\, \mathbf{h}^\top \mathbf{S}(\mathbf{a})}{\bar{g}^4}
+ \frac{12\, \mathbf{h}^\top \mathbf{I}\, \mathbf{h}}{\bar{g}^5} \right], \\
\frac{\partial^2 C_{ab}}{\partial a_k\, \partial h_0} &= \varepsilon \left[ \frac{\ell_k}{\bar{g}^2}
- \frac{2\, \mathbf{h}^\top \boldsymbol{\lambda}_k}{\bar{g}^3} \right], &
\frac{\partial^2 C_{ab}}{\partial h_0\, \partial \mathbf{h}} &= \varepsilon \left[ \frac{2\, \mathbf{S}(\mathbf{a})}{\bar{g}^3}
- \frac{6\, \mathbf{I}\, \mathbf{h}}{\bar{g}^4} \right], \\
\frac{\partial^2 C_{ab}}{\partial a_k\, \partial \mathbf{h}} &= \frac{\varepsilon\, \boldsymbol{\lambda}_k}{\bar{g}^2}, &
\frac{\partial^2 C_{ab}}{\partial \mathbf{h}\, \partial \mathbf{h}^\top} &= \frac{2 \varepsilon\, \mathbf{I}}{\bar{g}^3},
\end{aligned} \tag{57}
$$

and they follow in the degrees of freedom of the two bodies from (46), with
$\mathbf{J}_A = \mathbf{J}(\mathbf{c}_A)$ and
$\mathbf{J}_B = \mathbf{J}(\mathbf{c}_B)$:

$$
\frac{\partial C_{ab}}{\partial \mathbf{q}_A} = \mathbf{J}_A^\top \frac{\partial C_{ab}}{\partial \mathbf{s}}, \qquad
\frac{\partial C_{ab}}{\partial \mathbf{q}_B} = -\mathbf{J}_B^\top \frac{\partial C_{ab}}{\partial \mathbf{s}}, \qquad
\frac{\partial^2 C_{ab}}{\partial \mathbf{q}_A\, \partial \mathbf{q}_B^\top}
= -\mathbf{J}_A^\top \frac{\partial^2 C_{ab}}{\partial \mathbf{s}\, \partial \mathbf{s}^\top} \mathbf{J}_B, \tag{58}
$$

with the two diagonal blocks of the Hessian formed the same way and carrying
a positive sign. At rest, for body $A$ over a fixed handle and with moments
and edge integrals taken about $\mathbf{c}_A$, so that $d_x = d_y = 0$, the
gradient is

$$
\begin{aligned}
\frac{\partial C_{ab}}{\partial u_x} &= -\frac{\varepsilon\, \ell_x}{g_0}, &
\frac{\partial C_{ab}}{\partial u_y} &= -\frac{\varepsilon\, \ell_y}{g_0}, &
\frac{\partial C_{ab}}{\partial \theta_z} &= -\frac{\varepsilon\, \ell_\psi}{g_0}, \\
\frac{\partial C_{ab}}{\partial u_z} &= -\frac{\varepsilon A}{g_0^2}, &
\frac{\partial C_{ab}}{\partial \theta_x} &= -\frac{\varepsilon S_y}{g_0^2} + \frac{\varepsilon\, \ell_y d_z}{g_0}, &
\frac{\partial C_{ab}}{\partial \theta_y} &= \frac{\varepsilon S_x}{g_0^2} - \frac{\varepsilon\, \ell_x d_z}{g_0},
\end{aligned} \tag{59}
$$

where the terms in $d_z = z_0 - c_z$ arise because tilt about a center of
mass above or below the gap carries the face sideways. For a plate with no
handle outline under its silicon, $\Gamma_{ab} = \emptyset$, only gap and
tilt act, and the Hessian at rest is

$$
\frac{\partial^2 C_{ab}}{\partial (u_z, \theta_x, \theta_y)^2} = \frac{2 \varepsilon}{g_0^3}
\begin{bmatrix} A & S_y & -S_x \\ S_y & I_{yy} & -I_{xy} \\ -S_x & -I_{xy} & I_{xx} \end{bmatrix}. \tag{60}
$$

With potential difference $V$ between the two electrodes the co-energy of
the plate is $\tfrac{1}{2} C_{ab} V^2$. It applies generalized forces
$\tfrac{1}{2} V^2\, \partial C_{ab} / \partial \mathbf{q}$ to the two bodies
and adds $-\tfrac{1}{2} V^2\, \partial^2 C_{ab} / \partial \mathbf{q}\,
\partial \mathbf{q}^\top$ to their stiffness. Both come from the one function
(55), so force and stiffness agree with each other at every displacement.

!!! note "Implementation choices"
    The moments of $G_{ab}$ are found without building it. The sums of (11b)
    are evaluated for $U_{ab}$, for the oxide polygons inside it, and for the
    release openings inside it, and the last two are subtracted. Oxide and
    openings are disjoint, since (2) keeps oxide at least $d$ away from every
    opening, and each set is disjoint in itself. A polygon that $U_{ab}$
    covers is used as drawn and one that crosses its boundary is clipped to
    it first. A plate with a hundred thousand openings therefore costs one
    pass over their vertices, as for mass.

    $\Gamma_{ab}$ is built edge by edge. Each straight edge of each ring of
    the handle polygon is intersected with the silicon of the device
    electrode, and the openings that touch a resulting piece are cut out of
    it. A piece keeps the normal of the edge it came from, so no orientation
    has to survive the overlay: a ring edge from $(x_i, y_i)$ to
    $(x_{i+1}, y_{i+1})$ has the normal
    $(y_{i+1} - y_i,\; x_i - x_{i+1})$ to its right, which points out of the
    silicon for a counterclockwise exterior ring and for a clockwise hole,
    and is reversed otherwise. Orientation of a ring is the sign of its area
    sum in (11b).

    The region of a device electrode is the union of its body polygon and
    its stubs. Silicon of beams and bridges belongs to no electrode, so its
    capacitance to the handle is not counted.

    When both electrodes move with one body, (46) gives $\mathbf{s} = 0$ for
    every motion and the plate is a fixed capacitor. A bonded pad is the
    common case: oxide under its middle and the undercut rim of width $d$
    around it, both between silicon of the anchor.

    Fringing fields at the outline of the overlap and into release openings
    are not modeled. They add capacitance, of relative order $g_0 / w$ times
    a logarithmic factor for a plate of width $w$, and openings not much
    wider than the gap recover part of the area they remove, so (43) is a
    lower bound.

    The reference point is the center of the bounding box of the overlap
    rather than the origin or a center of mass, so that $\mathbf{S}$,
    $\mathbf{I}$, and the edge integrals stay of the size of the plate and
    the products in (58) do not cancel. In-plane travel is taken at
    mid-height of the gap for both bodies: taking each at its own face would
    turn a rigid tilt common to the two into travel of $t_\text{box}$ times
    that tilt. A rise that closes the gap at the reference point,
    $\bar{g} \le 0$, is an error.

## Sidewall capacitance

Facing sidewalls of device silicon form capacitors across trenches. Where the
device layer is thick against the width of a trench, the field between
sidewalls does not vary through the thickness, and capacitance is that of a
two-dimensional problem in the layout plane times the thickness. This section
treats the plane problem for a set of conductors: regions bounded by
polygons, conductor $a$ at potential $V_a$, in a medium of permittivity
$\varepsilon = \varepsilon_0 \varepsilon_a$. Every charge and capacitance
here is per unit thickness.

Charge sits on the outlines $\Gamma_a$ with density $\sigma$ and produces
the potential

$$
\varphi(\mathbf{r}) = -\frac{1}{2 \pi \varepsilon} \sum_a \oint_{\Gamma_a} \sigma(\mathbf{r}') \ln |\mathbf{r} - \mathbf{r}'|\, ds' + \varphi_0. \tag{61}
$$

In two dimensions the potential of a net charge grows as the logarithm of
distance, so capacitance to a ground at infinity does not exist. The set of
conductors is therefore neutral as a whole, and the constant $\varphi_0$,
the potential far away, is an unknown that floats:

$$
\varphi(\mathbf{r}) = V_a \quad \text{for } \mathbf{r} \text{ on } \Gamma_a, \qquad
\sum_a \oint_{\Gamma_a} \sigma\, ds = 0. \tag{62}
$$

Only differences of potential matter. Equations (61) and (62) have one
solution. Its potential is constant inside every conductor, so $\sigma$ is
the charge density that the field outside leaves on the surface.

### Elements

Outlines are divided into $N$ straight elements. Element $j$ has length
$L_j$ and carries charge $q_j$, spread evenly along it. Requiring the first
condition of (62) on average over each element, which is Galerkin's method
with constant weights, gives

$$
\sum_{j=1}^{N} P_{ij}\, q_j + \varphi_0 = V_{a(i)}, \quad i = 1, \dots, N, \qquad
\sum_{j=1}^{N} q_j = 0, \tag{63}
$$

where $a(i)$ is the conductor of element $i$ and

$$
P_{ij} = -\frac{I_{ij}}{2 \pi \varepsilon\, L_i L_j}, \qquad
I_{ij} = \int_{\text{element } i} \int_{\text{element } j} \ln |\mathbf{r} - \mathbf{r}'|\, ds'\, ds. \tag{64}
$$

$P_{ij}$ is the average potential on element $i$ of a unit charge on element
$j$. It is symmetric, which makes the capacitance matrix symmetric exactly
and not only in the limit of a fine mesh.

### Pair integral

The double integral of (64) has a closed form. Write a point as the complex
number $z = x + \mathrm{i} y$ and element $i$ as $z_i(s) = z_i^0 + \tau_i s$
for $0 \le s \le L_i$, with end points $z_i^0$ and $z_i^1$ and unit tangent
$\tau_i = (z_i^1 - z_i^0) / L_i$. Since $\ln |z| = \operatorname{Re} \ln z$,

$$
I_{ij} = \operatorname{Re} \int_0^{L_i} \int_0^{L_j} \ln \bigl( z_i(s) - z_j(s') \bigr)\, ds'\, ds. \tag{65}
$$

The argument of the logarithm is linear in $s$ and $s'$. With

$$
G(u) = u^2 \left( \tfrac{1}{2} \ln u - \tfrac{3}{4} \right), \qquad G''(u) = \ln u, \tag{66}
$$

$\partial^2 G(z_i(s) - z_j(s')) / \partial s\, \partial s' = -\tau_i \tau_j
\ln(z_i(s) - z_j(s'))$, so the integral is a second difference of $G$ over
the four pairs of end points, $u_{kl} = z_i^k - z_j^l$:

$$
I_{ij} = \operatorname{Re} \left\{ -\frac{1}{\tau_i \tau_j} \bigl[ G(u_{11}) - G(u_{10}) - G(u_{01}) + G(u_{00}) \bigr] \right\}. \tag{67}
$$

The four values $u_{kl}$ are the corners of a parallelogram in the complex
plane, and the logarithm must be one continuous branch over it. Two
different elements that do not cross leave the origin outside the
parallelogram or at one of its corners, so the branch cut can be laid away from it: every $u_{kl}$ is
multiplied by $\bar{u}_c / |u_c|$ before the principal logarithm is taken,
where $u_c$ is the difference of the element centers. The constant this
removes from the logarithm multiplies $L_i L_j$ by a purely imaginary number
and drops out of the real part. Elements that share an end point have
$u_{kl} = 0$ at one corner, where $G(0) = 0$. For an element with itself
(67) reduces to

$$
I_{ii} = L_i^2 \left( \ln L_i - \tfrac{3}{2} \right). \tag{68}
$$

Equation (67) subtracts numbers of size $|u|^2 \ln |u|$ to leave one of size
$L_i L_j \ln |u|$, so it loses accuracy when the elements are short against
their distance. For such pairs the logarithm is expanded about $u_c$
instead. With $\alpha = \tau_i L_i / (2 u_c)$ and
$\beta = \tau_j L_j / (2 u_c)$, odd powers integrate to zero and

$$
I_{ij} = L_i L_j \operatorname{Re} \left[ \ln u_c - \sum_{m=1}^{M} \frac{S_m}{2m} \right], \qquad
S_m = \sum_{k=0}^{m} \binom{2m}{2k} \frac{\alpha^{2k} \beta^{2m-2k}}{(2k+1)(2m-2k+1)}. \tag{69}
$$

The series is used when $L_i + L_j \le 2 \rho\, |u_c|$ with $\rho = 1/8$,
and is cut at $M = 7$. Its remainder is below $\rho^{16} / 16$, about
$2 \times 10^{-16}$ of $L_i L_j$.

### Capacitance matrix

Let $\mathbf{E}$ be the $N \times n$ matrix with $E_{ia} = 1$ when element
$i$ lies on conductor $a$ and zero otherwise, $\mathbf{v}$ the potentials of
the $n$ conductors, and $\mathbf{1}$ the vector of ones. Eliminating
$\varphi_0$ from (63) gives the element charges and the capacitance matrix,

$$
\mathbf{q} = \mathbf{M} \mathbf{E} \mathbf{v}, \qquad
\mathbf{M} = \mathbf{P}^{-1} - \frac{\mathbf{P}^{-1} \mathbf{1} \mathbf{1}^\top \mathbf{P}^{-1}}{\mathbf{1}^\top \mathbf{P}^{-1} \mathbf{1}}, \qquad
\mathbf{c} = \mathbf{E}^\top \mathbf{M} \mathbf{E}, \tag{70}
$$

where $c_{ab}$ is the charge on conductor $a$ per volt on conductor $b$ with
all others at zero. It has three properties by construction:

$$
\mathbf{c}^\top = \mathbf{c}, \qquad \mathbf{c}\, \mathbf{1} = \mathbf{0}, \qquad
\mathbf{v}^\top \mathbf{c}\, \mathbf{v} \ge 0. \tag{71}
$$

The second states that a common potential on all conductors moves no charge,
so $-c_{ab}$ is the capacitance between conductors $a$ and $b$. Stored
energy is $\tfrac{1}{2} \mathbf{v}^\top \mathbf{c}\, \mathbf{v}$. It equals
the largest value of $\mathbf{v}^\top \mathbf{E}^\top \mathbf{q} -
\tfrac{1}{2} \mathbf{q}^\top \mathbf{P} \mathbf{q}$ over neutral charges
$\mathbf{q}$, and charge that is constant on each element is a restricted
choice, so the computed energy never exceeds the exact one for the same
polygons and rises toward it as the mesh is refined. Capacitance of a layer
of thickness $t$ is $t\, \mathbf{c}$.

### Mesh

Charge density is smooth along a straight sidewall, diverges at a convex
corner of a conductor, as $r^{-1/3}$ at a right angle, and vanishes at a
concave one. Elements are therefore at most $h$ long and shorten toward
convex corners. Let a vertex turn the outline by the angle $\theta$ toward
the conductor, so that $\theta > 0$ at a convex corner and $\theta = \pi/2$
at a right angle. The length asked for at the vertex is

$$
h_c = h\, g^{-w}, \qquad w = \min \bigl( \max( 2 \theta / \pi,\, 0 ),\, 1 \bigr), \tag{72}
$$

which is $h / g$ at a right angle or sharper, and $h$ at a straight or
concave vertex, so the many vertices of a drawn arc cause no refinement.
Along an edge of length $\ell$ between vertices that ask for $h_L$ and
$h_R$, the length asked for at distance $x$ from the first grows with
distance from either end,

$$
s(x) = \min ( h,\; h_L + x,\; h_R + \ell - x ). \tag{73}
$$

The edge gets $n$ elements with nodes $x_k$ at equal steps of the element
count that $s$ implies,

$$
\nu(x) = \int_0^x \frac{dx'}{s(x')}, \qquad n = \lceil \nu(\ell) \rceil, \qquad
\nu(x_k) = \frac{k}{n}\, \nu(\ell), \quad k = 0, \dots, n. \tag{74}
$$

On a rising ramp $\nu = \ln(1 + x / h_L)$, so lengths there grow by a factor
of at most $e$ from one element to the next, and no element is longer than
$h$. An edge shorter than the lengths asked for at its ends is one element.

### Rigid motion

Conductors are grouped into rigid bodies. Body $B$ travels by $(a_x, a_y)$
and turns by $a_\psi$ about a point $\mathbf{c}_B$, written $c_B$ as a
complex number, which to first order moves a point $z$ of its outlines by

$$
\delta z = a_x + \mathrm{i}\, a_y + \mathrm{i}\, a_\psi\, (z - c_B). \tag{75}
$$

Elements keep their lengths, so $\mathbf{E}$ and $\mathbf{1}$ in (70) do not
change and the capacitance matrix depends on the motion through
$\mathbf{P}$ alone. $\mathbf{M}$ is the leading block of the inverse of the
bordered matrix of (63), whose border is constant, so for any parameter $p$
of the motion

$$
\frac{\partial \mathbf{c}}{\partial p} = -\mathbf{X}^\top \frac{\partial \mathbf{P}}{\partial p}\, \mathbf{X}, \qquad
\mathbf{X} = \mathbf{M} \mathbf{E}, \tag{76}
$$

where column $a$ of $\mathbf{X}$ holds the element charges for one volt on
conductor $a$. No further system is solved. The pair integral (65) depends
only on how two elements lie against each other, so its derivative is zero
when both belong to one body. For element $i$ on the moving body and
element $j$ on another, differentiating (65) under the integral with (75)
brings in

$$
T_{ij} = \int_0^{L_i} \int_0^{L_j} \frac{ds'\, ds}{z_i(s) - z_j(s')}, \qquad
W_{ij} = \int_0^{L_i} \int_0^{L_j} \frac{z_i(s) - c_B}{z_i(s) - z_j(s')}\, ds'\, ds, \tag{77}
$$

and gives

$$
\frac{\partial I_{ij}}{\partial a_x} = \operatorname{Re} T_{ij}, \qquad
\frac{\partial I_{ij}}{\partial a_y} = -\operatorname{Im} T_{ij}, \qquad
\frac{\partial I_{ij}}{\partial a_\psi} = -\operatorname{Im} W_{ij}. \tag{78}
$$

Both integrals are second differences like (67). With
$H(u) = G'(u) = u\, (\ln u - 1)$, the same branch of the logarithm, signs
$\epsilon_{kl} = +1$ for $k = l$ and $-1$ otherwise, and $I^\mathrm{c}_{ij}$
the complex number in braces in (67),

$$
T_{ij} = -\frac{1}{\tau_i \tau_j} \sum_{k, l} \epsilon_{kl}\, H(u_{kl}), \qquad
W_{ij} = -\frac{1}{\tau_i \tau_j} \sum_{k, l} \epsilon_{kl}\, H(u_{kl})\, (z_i^k - c_B) - I^\mathrm{c}_{ij}. \tag{79}
$$

For separated pairs the expansion that gave (69) gives, with $z_i^c$ the
center of element $i$,

$$
T_{ij} = \frac{L_i L_j}{u_c} \left[ 1 + \sum_{m=1}^{M} S_m \right], \qquad
W_{ij} = (z_i^c - c_B)\, T_{ij} - L_i L_j \sum_{m=1}^{M} R_m, \qquad
R_m = \sum_{k=1}^{m} \binom{2m-1}{2k-1} \frac{\alpha^{2k} \beta^{2m-2k}}{(2k+1)(2m-2k+1)}, \tag{80}
$$

in which the first term of $W_{ij}$ is the pull on the element acting at
its center and the second its turning about that center. Since
$I_{ij} = I_{ji}$, a pair whose second element moves is covered by swapping
the two. With $s_{ia} = X_{ia} / L_i$ the charge density on element $i$ per
volt on conductor $a$, (64) and (76) combine to

$$
\frac{\partial c_{ab}}{\partial p_B} = \frac{1}{2 \pi \varepsilon} \sum_{i \in B} \sum_{j \notin B}
\left( s_{ia}\, s_{jb} + s_{ib}\, s_{ja} \right) \frac{\partial I_{ij}}{\partial p_B} \tag{81}
$$

for each of the three parameters $p_B$ of body $B$.

Equation (81) is the derivative of the computed matrix (70) itself, not a
separate approximation of the derivative of the exact one, so force and
energy belong to one function of position: with potentials $\mathbf{v}$ the
force on body $B$ is $\tfrac{1}{2} \mathbf{v}^\top (\partial \mathbf{c} /
\partial p_B)\, \mathbf{v}$ per unit thickness, and it equals the change of
$\tfrac{1}{2} \mathbf{v}^\top \mathbf{c}\, \mathbf{v}$ when the bodies and
their elements are moved, to rounding. It is the pull of the charge on other
bodies on the charge of body $B$. The field of other bodies is smooth along
the outlines of $B$, so the squares of charge density that diverge at
corners, which a surface stress integral would need, never appear. Moving
every body together changes nothing:

$$
\sum_B \frac{\partial \mathbf{c}}{\partial a_x^B} = \sum_B \frac{\partial \mathbf{c}}{\partial a_y^B} = \mathbf{0}, \qquad
\sum_B \frac{\partial \mathbf{c}}{\partial a_\psi^B} = \mathbf{0} \quad \text{for a common } \mathbf{c}_B, \tag{82}
$$

which holds pair by pair, because $T_{ij} + T_{ji} = 0$ and
$W_{ij} + W_{ji} = L_i L_j$ is real: forces and torques on the bodies
balance.

!!! note "Implementation choices"
    Lengths are divided by twice the diagonal of the bounding box of all
    outlines before (64) is formed. Every distance is then below one, which
    makes $\mathbf{P}$ positive definite, and (70) is solved with one
    Cholesky factorization and $n + 1$ right-hand sides. The scaling adds
    the same constant to every $P_{ij}$, which cancels through the second
    condition of (63), and $\varepsilon$ is applied to the result.

    A hole in a conductor is meshed only when another polygon lies inside
    it. A cavity that holds no conductor holds no field and its wall carries
    no charge.

    Conductors must not overlap or touch: two regions in contact at
    different potentials have no finite capacitance. Touching regions at
    one potential are passed as one conductor.

    In (67) the remaining loss of accuracy is set by how short an element is
    against a neighbor, which (72) bounds by $g$. Measured against 40-digit
    arithmetic over random pairs with length ratios up to 64, both forms
    agree with the exact integral to $10^{-12}\, L_i L_j$.

    The system is dense: memory grows as $N^2$ and time as $N^3$. A mesh
    with more elements than `max_elements` is an error rather than a solve
    that does not finish. The gradient adds one pass over the pairs of
    elements that lie on different bodies.

    A conductor couples to every other conductor in the plane, however far,
    because the plane field decays only as the logarithm of distance. Which
    conductors enter one solution is the caller's choice.

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
| $\rho_e$ | electrical resistivity of the device layer | `Material.resistivity` |
| $R_v$ | resistance of one via link | `ExtractionParams.via_resistance` |
| $g_0$ | gap between device and handle silicon | `Layer.thickness` of the box layer |
| $\varepsilon_\text{ox}$ | relative permittivity of buried oxide | `Material.relative_permittivity` of the box layer |
| $\varepsilon_a$ | relative permittivity of the medium in released gaps and trenches | `ElectrostaticParams.ambient_permittivity` |
| $h$ | longest boundary element | `ElectrostaticParams.element_size` |
| $g$ | refinement of elements at convex corners | `ElectrostaticParams.corner_grading` |
| $E$, $G$ | Young's and shear modulus of the device layer | `Material.youngs_modulus`, `Material.shear_modulus` or `Material.poisson_ratio` |
| $J$ | torsion constant model | `Material.torsion_model` |
| $\beta$ | short-beam torsional stiffening coefficient | `Material.torsion_beta` |
