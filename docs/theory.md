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
conductors joined through vias and through wires: connections made off the
die, such as wire bonds between pads, bonds to a package, or traces on a
board, which the layout does not show and the caller supplies.

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
- A wire joins the electrodes that cover its two ends. Each end is a point
  on a pad of a named layer, and the conductor that covers it is the one the
  wire joins in the nets. Its link has the resistance given with the wire,
  which stands for the bond or trace and its two contacts.

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
one conductor, and a via or a wire joins conductors of one net. Conversely, the
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

    A wire end must lie on silicon of an electrode. An end on a beam or a
    bridge is an error, since that silicon belongs to no electrode, and so
    is an end on no conductor of its layer, which also catches silicon that
    pruning dropped. A point in a release hole of a pad counts as on the
    pad, as the pad's region is filled. A wire with both ends on one
    electrode joins nothing and is an error.

    Wires belong to the nets, so every analysis that reads nets sees bonded
    pads as one: the capacitance between nets, potentials and charges in
    static equilibrium, and the resistor network of the links.

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

Two methods find the capacitance matrix of the plane problem. The field
method of this section solves the field between the conductors. It
approaches the exact answer of the plane problem as its mesh is refined,
and a mesh too large for one solution is solved in
[windows](#windows). The plates method of the
[next section](#sidewalls-as-facing-plates) takes every wall that faces
another conductor as a parallel plate with it. It costs far less and leaves
out all fringe field. Which one suits a
device is a judgment of its geometry, so `sidewall_method` has no default:
the caller chooses.

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

### Second derivatives

Electrostatic stiffness needs the second derivative of capacitance, for
which (75) must be the rigid motion it approximates. Body $B$ travels by
$a = a_x + \mathrm{i}\, a_y$ and turns by $a_\psi$ about a point that travels
with it:

$$
z \mapsto c_B + a + e^{\mathrm{i} a_\psi} (z - c_B). \tag{83}
$$

Differentiating (76) once more, with $\partial \mathbf{M} / \partial q =
-\mathbf{M}\, (\partial \mathbf{P} / \partial q)\, \mathbf{M}$, gives

$$
\frac{\partial^2 \mathbf{c}}{\partial p\, \partial q}
= -\mathbf{X}^\top \frac{\partial^2 \mathbf{P}}{\partial p\, \partial q} \mathbf{X}
+ \mathbf{Y}_p^\top \mathbf{M} \mathbf{Y}_q + \mathbf{Y}_q^\top \mathbf{M} \mathbf{Y}_p, \qquad
\mathbf{Y}_p = \frac{\partial \mathbf{P}}{\partial p} \mathbf{X}. \tag{84}
$$

The first term is the change of the pull between fixed charges. The other
two are the pull that arises because the charges themselves move along the
outlines when a body moves; they need $\mathbf{M} \mathbf{Y}_q$, one more
pair of triangular solves per parameter with the factor already at hand.

Take a pair with element $i$ on body $B$ and element $j$ on body $B'$, and
the six parameters $(a_x, a_y, a_\psi)$ of $B$ and $(a_x', a_y', a_\psi')$
of $B'$. With $r_k = z_i^k - c_B$ and $r_l' = z_j^l - c_{B'}$, the end-point
differences $u_{kl}$ of (67) move under (83) with

$$
\frac{\partial u_{kl}}{\partial p} = v_p = \bigl( 1,\; \mathrm{i},\; \mathrm{i}\, r_k,\; -1,\; -\mathrm{i},\; -\mathrm{i}\, r_l' \bigr), \qquad
\frac{\partial^2 u_{kl}}{\partial a_\psi^2} = -r_k, \qquad
\frac{\partial^2 u_{kl}}{\partial a_\psi'^2} = r_l', \tag{85}
$$

all other second derivatives being zero, and the product $\tau_i \tau_j$
turns by $e^{\mathrm{i} (a_\psi + a_\psi')}$. Let $\varrho_p$ be one for the
two rotations and zero for the four translations. Differentiating the
complex integral $I^\mathrm{c}_{ij}$ of (67) gives

$$
\begin{aligned}
\frac{\partial I^\mathrm{c}_{ij}}{\partial p} &= \Phi_p - \mathrm{i}\, \varrho_p\, I^\mathrm{c}_{ij}, \qquad
\Phi_p = -\frac{1}{\tau_i \tau_j} \sum_{k,l} \epsilon_{kl}\, H(u_{kl})\, v_p, \\
\frac{\partial^2 I^\mathrm{c}_{ij}}{\partial p\, \partial q} &= -\frac{1}{\tau_i \tau_j} \sum_{k,l} \epsilon_{kl}
\left[ \ln (u_{kl})\, v_p v_q + H(u_{kl}) \frac{\partial^2 u_{kl}}{\partial p\, \partial q} \right]
- \mathrm{i}\, \varrho_p \Phi_q - \mathrm{i}\, \varrho_q \Phi_p - \varrho_p \varrho_q\, I^\mathrm{c}_{ij},
\end{aligned} \tag{86}
$$

whose real parts are the derivatives of $I_{ij}$. The first line contains
(78) and (79): $\Phi_{a_x} = T_{ij}$ and
$\Phi_{a_\psi} - \mathrm{i} I^\mathrm{c}_{ij} = \mathrm{i}\, W_{ij}$.

For separated pairs write (69) as $I^\mathrm{c}_{ij} = L_i L_j\, F$ with
$F = \ln u_c - \sum_{m,k} e_{mk}\, t_{mk} / (2m)$, where
$t_{mk} = \alpha^{2k} \beta^{2m-2k}$ and $e_{mk}$ is its coefficient in
$S_m$. $F$ depends on the motion through $u_c$ and through the two angles,
since $\alpha^2$ turns by $e^{2 \mathrm{i} a_\psi}$ and $\beta^2$ by
$e^{2 \mathrm{i} a_\psi'}$. Its partial derivatives are sums over the same
terms,

$$
\begin{aligned}
F_u &= \frac{1}{u_c} \Bigl[ 1 + \sum e_{mk} t_{mk} \Bigr], &
F_{uu} &= -\frac{1}{u_c^2} \Bigl[ 1 + \sum (2m+1)\, e_{mk} t_{mk} \Bigr], \\
F_\psi &= -\mathrm{i} \sum \frac{k}{m}\, e_{mk} t_{mk}, &
F_{\psi'} &= -\mathrm{i} \sum \frac{m-k}{m}\, e_{mk} t_{mk}, \\
F_{u\psi} &= \frac{2 \mathrm{i}}{u_c} \sum k\, e_{mk} t_{mk}, &
F_{u\psi'} &= \frac{2 \mathrm{i}}{u_c} \sum (m-k)\, e_{mk} t_{mk}, \\
F_{\psi\psi} &= \sum \frac{2 k^2}{m}\, e_{mk} t_{mk}, &
F_{\psi'\psi'} &= \sum \frac{2 (m-k)^2}{m}\, e_{mk} t_{mk}, \qquad
F_{\psi\psi'} = \sum \frac{2 k (m-k)}{m}\, e_{mk} t_{mk}.
\end{aligned} \tag{87}
$$

The difference of centers $u_c$ moves as in (85) with the element centers in
place of the end points; call its first derivatives $w_p$. Writing
$\omega_p$ for the angle that parameter $p$ is, if it is one,

$$
\frac{\partial F}{\partial p} = F_u\, w_p + F_{\omega_p}, \qquad
\frac{\partial^2 F}{\partial p\, \partial q} = F_{uu}\, w_p w_q + F_u \frac{\partial^2 u_c}{\partial p\, \partial q}
+ F_{u \omega_p} w_q + F_{u \omega_q} w_p + F_{\omega_p \omega_q}, \tag{88}
$$

where a term with $\omega_p$ is absent when $p$ is a translation. The first
of these contains (80).

With these, the first term of (84) is a sum over pairs of elements on
different bodies, as (81) is,

$$
-\Bigl[ \mathbf{X}^\top \frac{\partial^2 \mathbf{P}}{\partial p\, \partial q} \mathbf{X} \Bigr]_{ab}
= \frac{1}{2 \pi \varepsilon} \sum_{B < B'} \sum_{i \in B} \sum_{j \in B'}
\left( s_{ia}\, s_{jb} + s_{ib}\, s_{ja} \right) \frac{\partial^2 I_{ij}}{\partial p\, \partial q}, \tag{89}
$$

in which a pair contributes when $p$ and $q$ are each a parameter of one of
its two bodies, and

$$
\bigl[ \mathbf{Y}_p \bigr]_{ia} = -\frac{1}{2 \pi \varepsilon\, L_i} \sum_j \frac{\partial I_{ij}}{\partial p}\, s_{ja}, \tag{90}
$$

summed over the elements $j$ that do not share a body with element $i$.
Translating every body together leaves the gradient unchanged, so the
second derivatives for travel along one axis sum to zero over the bodies of
either parameter.

### Windows

One dense solution of (70) is out of reach for a comb of thousands of
fingers. Charge on a sidewall, though, is set by its surroundings: across a
narrow gap by the wall it faces, and hardly at all by conductors far off. A
mesh with more elements than `max_elements` is therefore solved piece by
piece.

The elements are divided into cores, every element in exactly one, by
halving the set at the median of the element centers along the longer side
of their bounding box until no core holds more than $n_w$ elements. Let
$K_k$ be the elements of core $k$ and $B_k$ the bounding box of their
centers $\mathbf{r}_i$. Window $k$ holds the core and every element within
the overlap $\delta_k$ of that box,

$$
W_k = \bigl\{ i : \mathbf{r}_i \in B_k \oplus [-\delta_k, \delta_k]^2 \bigr\} \supseteq K_k, \tag{91}
$$

and is solved as a problem (63) of its own: its elements alone, each at the
potential of its conductor, neutral together, with a floating potential of
its own. Let $\mathbf{E}_k$ and $\mathbf{M}_k$ be the matrices of (70) for
the elements of $W_k$, and $\mathbf{E}_k^c$ equal to $\mathbf{E}_k$ with the
rows of elements outside the core set to zero. The window supplies the
charge on its core alone,

$$
\mathbf{c}^{(k)} = \mathbf{E}_k^{c\,\top} \mathbf{M}_k \mathbf{E}_k, \tag{92}
$$

whose entry $(a, b)$ is the charge on the core elements of conductor $a$ per
volt on conductor $b$. Every element is in one core, so the sum over windows
counts the charge on every element once:

$$
\hat{\mathbf{c}} = \sum_k \mathbf{c}^{(k)}, \qquad \hat{\mathbf{c}}\, \mathbf{1} = \mathbf{0}. \tag{93}
$$

Rows sum to zero because $\mathbf{M}_k \mathbf{1} = \mathbf{0}$ in every
window. $\hat{\mathbf{c}}$ is not symmetric, since each window sees its own
part of the field, and the capacitance matrix is

$$
c_{ab} = \frac{\hat{c}_{ab} + \hat{c}_{ba}}{2} \quad (a \ne b), \qquad
c_{aa} = -\sum_{b \ne a} c_{ab}, \tag{94}
$$

which has the first two properties of (71). The asymmetry
$\max |\hat{c}_{ab} - \hat{c}_{ba}| / \max |\hat{c}_{ab}|$ is zero when one
window holds the whole mesh. Otherwise it measures what the windows leave
out, and it is reported with the result.

Elements keep their cores and windows as bodies move, so (92) depends on
the motion through $\mathbf{P}_k$ alone, as (70) does. With
$\mathbf{X}_k = \mathbf{M}_k \mathbf{E}_k$ and
$\mathbf{Y}_k = \mathbf{M}_k \mathbf{E}_k^c$,

$$
\begin{aligned}
\frac{\partial \mathbf{c}^{(k)}}{\partial p} &= -\mathbf{Y}_k^\top \frac{\partial \mathbf{P}_k}{\partial p}\, \mathbf{X}_k, \\
\frac{\partial^2 \mathbf{c}^{(k)}}{\partial p\, \partial q} &= -\mathbf{Y}_k^\top \frac{\partial^2 \mathbf{P}_k}{\partial p\, \partial q}\, \mathbf{X}_k
+ \Bigl( \frac{\partial \mathbf{P}_k}{\partial p} \mathbf{Y}_k \Bigr)^{\!\top} \mathbf{M}_k\, \frac{\partial \mathbf{P}_k}{\partial q} \mathbf{X}_k
+ \Bigl( \frac{\partial \mathbf{P}_k}{\partial q} \mathbf{Y}_k \Bigr)^{\!\top} \mathbf{M}_k\, \frac{\partial \mathbf{P}_k}{\partial p} \mathbf{X}_k,
\end{aligned} \tag{95}
$$

which are (76) and (84) with $\mathbf{Y}_k$ in place of $\mathbf{X}_k$ on
the left, evaluated with the same pair derivatives (85) to (88). They pass
through (93) and (94) as $\mathbf{c}^{(k)}$ does and are the exact
derivatives of the matrix (94).

The overlap is found for each window. It starts at $\delta_k = 4 h$ and
doubles until the charge on the core has settled,

$$
\max_{a,b} \bigl| c^{(k)}_{ab}(\delta_k) - c^{(k)}_{ab}(\delta_k / 2) \bigr| < \tau\, \max_{a,b} \bigl| c^{(k)}_{ab}(\delta_k) \bigr|, \tag{96}
$$

with $\tau$ the window tolerance, or until $\delta_k$ reaches the interaction
distance $d_i$ when one is given, or until the next window would hold more
than `max_elements`. A window whose elements all lie on one conductor
supplies nothing and is not solved; its overlap still doubles, since a
second conductor may come into view.

Measured against one dense solution of the same mesh, a comb of 30 fingers
with gaps of 3 um, 4825 elements in 32 windows:

| $d_i$ | Capacitance | Force along fingers | Stiffness across fingers | Asymmetry |
| --- | --- | --- | --- | --- |
| 30 um | $-2.6 \times 10^{-4}$ | $-6.6 \times 10^{-3}$ | $-4.1 \times 10^{-5}$ | $4.1 \times 10^{-4}$ |
| 60 um | $+3.2 \times 10^{-4}$ | $-7.4 \times 10^{-4}$ | $+2.3 \times 10^{-6}$ | $5.2 \times 10^{-5}$ |

Charge in the gaps settles within an overlap of two gaps. Force along the
fingers converges slowest: it comes from the charge at finger tips, which
face open space, where the plane field reaches far. $d_i$ belongs at a
few times the widest trench that should count: a window whose overlap
stops short of two trench widths leaves out part of the field across it.

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

    One solution is dense: memory grows as $N^2$ and time as $N^3$, and
    gradient and second derivatives add one pass over the pairs of elements
    that lie on different bodies, about as long as forming $\mathbf{P}$. In
    windows the cost grows with $N$ alone and memory with the size of one
    window, so a cluster far too large to hold as one system is solved in
    the time of its windows. Every window is scaled by its own bounding box.

    A window is neutral and floats, as the whole problem does, because no
    other condition at its cut edge is right in the plane: charge in a gap
    depends only on the difference of the potentials across it, which a
    floating window keeps, while charge on a wall that faces open space
    depends on the potential far away, which a window can only guess. That
    is the part $d_i$ and $\tau$ trade against time. When $d_i$ is given the
    overlap stops there, the distance beyond which the plane field is not
    trusted anyway.

    A first window that already holds more than `max_elements` is an error.

    The point a body turns about travels with the body, as a center of mass
    does. A second derivative for travel and rotation taken about a point
    fixed in space differs from this one by the gradient.

    A conductor couples to every other conductor in the plane, however far,
    because the plane field decays only as the logarithm of distance. Which
    conductors enter one solution is the caller's choice.

## Sidewalls as facing plates

In a comb most of the capacitance sits in narrow gaps between long facing
walls, where the field is that of a parallel plate. The plates method keeps
that part alone. It solves no field: every piece of wall that looks straight
across a gap at another conductor is a plate with it.

From a point $\mathbf{r}$ on the outline of conductor $a$, follow the
outward normal to the first outline it meets. Let $\gamma(\mathbf{r})$ be
the distance to that outline and $b(\mathbf{r})$ its conductor. With $d_i$
the interaction distance, capacitance per unit thickness between conductors
$a$ and $b$ is

$$
-c_{ab} = \frac{\varepsilon}{2} \left[ \int_{\Gamma_{a \to b}} \frac{ds}{\gamma} + \int_{\Gamma_{b \to a}} \frac{ds}{\gamma} \right], \qquad
\Gamma_{a \to b} = \bigl\{ \mathbf{r} \in \Gamma_a : b(\mathbf{r}) = b,\ \gamma(\mathbf{r}) < d_i \bigr\}, \tag{97}
$$

for $a \ne b$. Each integral lays plates of width $ds$ and gap $\gamma$
along the walls of one conductor that face the other. The two views of a gap
agree when its walls are parallel and differ in second order of the angle
between them otherwise, and their mean makes $c_{ab}$ symmetric. A wall
whose normal first meets its own conductor, as inside a slot or an etch
hole, contributes nothing, and neither does a wall that meets nothing within
$d_i$. The diagonal follows from

$$
c_{aa} = -\sum_{b \ne a} c_{ab}, \tag{98}
$$

so the matrix has the first two properties of (71), and $-c_{ab}$ is the
capacitance between two conductors as there.

### Pieces of wall

Outlines are polygons. Edge $i$ runs from $\mathbf{r}_i^0$ over the length
$L_i$ with unit tangent $\boldsymbol{\tau}_i$ and outward unit normal
$\mathbf{n}_i$, and a point $\mathbf{r}$ has the coordinates

$$
u = (\mathbf{r} - \mathbf{r}_i^0) \cdot \boldsymbol{\tau}_i, \qquad
v = (\mathbf{r} - \mathbf{r}_i^0) \cdot \mathbf{n}_i \tag{99}
$$

in the frame of the edge: $u$ along the wall and $v$ out of the silicon. The
normal from the point at $u$ is the line of constant $u$. Another edge $j$
with end points $(u_j^0, v_j^0)$ and $(u_j^1, v_j^1)$ in this frame meets
that line at the height

$$
v_j(u) = v_j^0 + m_j\, (u - u_j^0), \qquad m_j = \frac{v_j^1 - v_j^0}{u_j^1 - u_j^0}, \tag{100}
$$

for $u$ between $u_j^0$ and $u_j^1$. The gap at $u$ is the lowest height
above the wall, and the edge that attains it is the edge faced:

$$
\gamma(u) = \min_j \bigl\{ v_j(u) : v_j(u) > 0 \bigr\}. \tag{101}
$$

Outlines do not cross, so the heights of two edges never cross between
their end points, and the edge faced can change only at a mark: the
projection $u_j^0$ or $u_j^1$ of an end point of some edge, a point where
$v_j(u) = d_i$, or an end $0$ or $L_i$ of the wall itself. The marks cut the
wall into intervals. On each, one edge is faced throughout or none lies
within $d_i$, and (101) at the middle of the interval tells which. Neighboring
intervals that face the same edge join into one piece $u_1 \le u \le u_2$,
and pieces that face the conductor of edge $i$ itself are dropped. Nothing
is sampled: a piece ends exactly where the geometry says, so capacitance is
a smooth function of position wherever no two marks coincide.

An edge that rises more than ten times its run along the wall,
$|m_j| > 10$, is no plate with the wall. Where it is the lowest edge it
hides what lies behind it, and that stretch of the wall counts for nothing.

### Integral over a piece

On a piece the gap is the linear function (100) of one edge $j$. With the run
$\Delta = u_2 - u_1$ and the gap at the middle
$\bar{\gamma} = v_j \bigl( (u_1 + u_2) / 2 \bigr)$,

$$
\int_{u_1}^{u_2} \frac{du}{v_j(u)} = \frac{\Delta}{\bar{\gamma}}\, \frac{\operatorname{artanh} x}{x}, \qquad
x = \frac{m_j\, \Delta}{2 \bar{\gamma}}. \tag{102}
$$

The first factor is the parallel plate of the mean gap. The second is
$1 + x^2 / 3 + x^4 / 5 + \dots$: one for parallel walls, and larger the
more one end of the piece closes. Walls that do not touch have $|x| < 1$.
The integrals of (97) are the sums of (102) over the pieces of all edges of
one conductor that face the other.

### Derivatives for rigid motion

A piece depends on the end points of at most four edges: edge $i$, which
sets the frame (99); the edge $j$ it faces, which sets (100); and the edges
whose end points or whose crossing of $d_i$ set $u_1$ and $u_2$. Each end
point moves with its body by (83), so (102) is a closed form of the
parameters $(a_x, a_y, a_\psi)$ of up to four bodies. Its derivatives come
from carrying, with every intermediate value $f$, the gradient $\nabla f$
and Hessian $\nabla \nabla f$ in those parameters through the arithmetic of
(99) to (102). A product and a function $\phi$ of one value follow

$$
\begin{aligned}
\nabla (f_1 f_2) &= f_1 \nabla f_2 + f_2 \nabla f_1, &
\nabla \nabla (f_1 f_2) &= f_1 \nabla \nabla f_2 + f_2 \nabla \nabla f_1 + \nabla f_1 \nabla f_2^\top + \nabla f_2 \nabla f_1^\top, \\
\nabla \phi(f) &= \phi'(f)\, \nabla f, &
\nabla \nabla \phi(f) &= \phi'(f)\, \nabla \nabla f + \phi''(f)\, \nabla f\, \nabla f^\top,
\end{aligned} \tag{103}
$$

with $\phi$ the reciprocal in a quotient and the second factor in (102). The
arithmetic starts at the end points, whose coordinates have by (83) the
gradient $(1, \mathrm{i}, \mathrm{i}\, (z - c_B))$ and the one second
derivative $-(z - c_B)$ for rotation twice. The result is the exact gradient
and Hessian of (97), in the same layout as (81) and (84) of the field
method.

When a third conductor stands between two walls, its corners set where the
pieces of those walls end. Capacitance between $a$ and $b$ then changes
with the motion of a body that carries neither, which (103) includes.

### What plates leave out

The plates method has no fringe field: none around finger tips, outer walls,
and corners, and none from a wall to anything but the first conductor it
faces. Measured against the field method with elements of 1 um, as the ratio
of plates to field with $d_i$ of 30 um:

| Geometry | Capacitance | Force | Stiffness |
| --- | --- | --- | --- |
| Comb, 8 gaps of 3 um, overlap 40 um | 0.86 | 1.01 along fingers | 0.96 across fingers |
| Wall of 200 um at a gap of 4 um | 0.89 | 0.98 | 0.99 |
| Racks of teeth 4 um wide on a pitch of 8 um, tops 2 um apart | 0.72 aligned, 0.36 at half a pitch | 4.3 at a quarter pitch | none, where the field has a detent |

Capacitance is low by the share of the fringe field. Force and stiffness of
a comb or a plate are close, because their fringe field hardly changes with
travel. A toothed motor is the opposite case: its force is the change of
the fringe field as teeth pass each other, and plates see only the overlap
of tooth tops. Such a device needs the field method.

Teeth that stand tip to tip at rest are the worst case. By plates their
capacitance is a cusp in travel: it falls to either side of rest, the two
one-sided forces are large and opposite, and the force at rest is whichever
of them the rounding of coordinates selects, and it can exceed the force
of the field many times over. The electrostatic spring there is the mean of
the two sides and is far less affected.

!!! note "Implementation choices"
    Edges in view of a wall are found in a spatial index of all edges,
    queried with the rectangle that the wall sweeps along its normal over
    $d_i$. An edge that lies along the normal, with end points less than
    $10^{-9} L_i$ apart in $u$, has no height (100) and is skipped. Marks
    closer than $10^{-9} L_i$ count as one.

    Below $|x| = 1/4$ the second factor of (102) is summed from its series,
    so nothing cancels as $m_j \to 0$.

    Where a wall lies exactly along the normal of another, as the end face
    of a finger does beside the wall its flank faces, (97) has a kink in
    rotation: a turn one way brings the end face into view and a turn the
    other way does not. Derivatives are those with it out of view. The step
    in $\partial c / \partial a_\psi$ is
    $\tfrac{1}{2} \varepsilon \ln (1 + w / \gamma)$ for an end face of width
    $w$ beside a gap $\gamma$, small against the
    $\varepsilon L r / \gamma^2$ of a flank of length $L$ at distance $r$
    from the pivot. Where marks set by end points of different bodies
    coincide, as when a finger tip lines up exactly with a corner across
    the gap, there is a kink in travel too. Of coincident marks the one that
    sets an end of a piece is an end point of the edge faced if there is
    one, else an end of the wall itself, else an end point of the edge
    faced on the other side of the mark.

    The limit on slope makes a wall that is square to another only up to
    the rounding of its coordinates behave as one that is square exactly.
    Taken as a plate it would be a sliver a few picometers wide, whose
    second derivative grows with its slope without bound, so on teeth turned
    slightly out of the grid such slivers would give an electrostatic spring
    of either sign and of any size.

    The cutoff is sharp. A wall parallel to another just inside $d_i$
    counts in full and just outside not at all, so $d_i$ belongs well above
    every gap that matters and below the distance to silicon that should
    not count.

    The method finds no charge density and needs no mesh. The pieces of
    wall, with the conductor each faces and the gap at both ends, are kept
    for review.

## Capacitance of a structure

The electrodes of a structure, $n$ in all, are the nodes of one capacitance
matrix $\mathbf{C}$, defined by the charges $\mathbf{Q} = \mathbf{C}
\mathbf{v}$ that potentials $\mathbf{v}$ put on them. It is the sum of the
plates between device and handle electrodes and of the sidewalls between
device electrodes.

A plate of capacitance $c$ between electrodes $a$ and $b$ stores energy
$\tfrac{1}{2} c\, (v_a - v_b)^2$, so with $\mathbf{e}_a$ the unit vector of
electrode $a$ it adds

$$
c\, (\mathbf{e}_a - \mathbf{e}_b)(\mathbf{e}_a - \mathbf{e}_b)^\top \tag{104}
$$

to $\mathbf{C}$: $c$ to the two diagonal entries and $-c$ to the two
others. Each plate enters twice, with $C_\text{ox}$ of (43) into
$\mathbf{C}^\text{ox}$ and with $\varepsilon A / g_0$ into
$\mathbf{C}^\text{rel}$.

Sidewalls are found in the plane of the device layer, by the method the
caller chooses. In that plane the field couples every conductor to every
other, however far, which the handle beneath and the space above prevent in
the real stack. Electrodes are therefore taken together only when their
silicon is close. Let $d_i$ be the interaction distance and join device
electrodes $a$ and $b$ when

$$
\operatorname{dist}(E_a, E_b) \le d_i. \tag{105}
$$

A cluster is a connected set of this relation with two or more electrodes:
electrodes joined directly or through a chain of neighbors. Each cluster is
one plane problem with its electrodes as conductors. With the field method
it is one problem (62), all conductors coupled to one another whatever
their distance, and $\mathbf{c}$ is that of (70), or that of (94) when the
cluster needs more elements than `max_elements` and is solved in windows
whose overlap stops at the same $d_i$. With the plates method $\mathbf{c}$
is that of (97) and (98), where $d_i$ is the widest gap across which two
walls count. Its matrix times the thickness of the device layer fills the
rows and columns of its electrodes,

$$
\bigl[ \mathbf{C}^\text{sw} \bigr]_{ab} = t_\text{device}\, c_{ab} \quad \text{for } a, b \text{ in one cluster}, \tag{106}
$$

and electrodes of different clusters share no sidewall capacitance. Then

$$
\mathbf{C} = \mathbf{C}^\text{ox} + \mathbf{C}^\text{rel} + \mathbf{C}^\text{sw}, \qquad
\mathbf{C}^\top = \mathbf{C}, \qquad \mathbf{C}\, \mathbf{1} = \mathbf{0}. \tag{107}
$$

Rows sum to zero because (104), (71), and (98) do: a potential common to all
electrodes moves no charge, nothing outside the structure holds charge, and
$-C_{ab}$ is the capacitance between electrodes $a$ and $b$.

Electrodes of one net are joined by links and share a potential at rest.
With $S_{a\nu} = 1$ when electrode $a$ belongs to net $\nu$ and zero
otherwise, the capacitance matrix between nets is

$$
\mathbf{C}^\text{net} = \mathbf{S}^\top \mathbf{C}\, \mathbf{S}, \tag{108}
$$

which adds the entries of $\mathbf{C}$ over the electrodes of each pair of
nets. Capacitance between two electrodes of one net drops out.

### Forces and stiffness

$\mathbf{C}$ depends on the positions of the bodies. Let $\mathbf{q}$ stack
the degrees of freedom (12) of every movable body, in the order of the
assembled matrices (34) and (35). A plate between two bodies contributes
its gradient and Hessian (58) in the pattern of (104). A cluster contributes
through the in-plane motion of its bodies at the mid-plane of the device
layer, $z_m = -t_\text{device} / 2$, where sidewall forces act. Body $B$
with center of mass $\mathbf{c}_B$ turns about the vertical through
$\mathbf{c}_B$ and, by (44), moves there by

$$
\begin{bmatrix} a_x \\ a_y \\ a_\psi \end{bmatrix} = \mathbf{J}_B\, \mathbf{q}_B, \qquad
\mathbf{J}_B = \begin{bmatrix}
1 & 0 & 0 & 0 & d_z & 0 \\
0 & 1 & 0 & -d_z & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1
\end{bmatrix}, \qquad d_z = z_m - c_{B,z}, \tag{109}
$$

so that with the gradient and second derivative of the cluster, (81) and
(84) for the field method or those of (103) for plates, taken about the
centers of mass,

$$
\frac{\partial \mathbf{C}^\text{sw}}{\partial \mathbf{q}_B} = t_\text{device}\, \frac{\partial \mathbf{c}}{\partial \mathbf{a}_B}\, \mathbf{J}_B, \qquad
\frac{\partial^2 \mathbf{C}^\text{sw}}{\partial \mathbf{q}_B\, \partial \mathbf{q}_{B'}^\top}
= t_\text{device}\, \mathbf{J}_B^\top \frac{\partial^2 \mathbf{c}}{\partial \mathbf{a}_B\, \partial \mathbf{a}_{B'}^\top}\, \mathbf{J}_{B'}. \tag{110}
$$

The lever $d_z$ is zero for a body of device silicon alone and large for
one that carries handle silicon: an in-plane pull on its sidewalls then
also tilts it. The anchor has no degrees of freedom and its rows and
columns are dropped, as in (37).

With potentials $\mathbf{v}$ on the electrodes the co-energy is
$\tfrac{1}{2} \mathbf{v}^\top \mathbf{C}(\mathbf{q})\, \mathbf{v}$, and the
generalized forces and the electrostatic stiffness are

$$
f_k = \frac{1}{2}\, \mathbf{v}^\top \frac{\partial \mathbf{C}}{\partial q_k}\, \mathbf{v}, \qquad
\bigl[ \mathbf{K}_e \bigr]_{kl} = -\frac{1}{2}\, \mathbf{v}^\top \frac{\partial^2 \mathbf{C}}{\partial q_k\, \partial q_l}\, \mathbf{v}. \tag{111}
$$

$\mathbf{K}_e$ adds to $\mathbf{K}^\text{full}$ of (35). It is symmetric,
and a negative diagonal entry means that the field softens that motion. For
a comb of $N$ gaps of width $g$ and overlap $\ell$ between a stator at $V$
and a shuttle at zero, the ideal values that (111) approaches are

$$
f_\parallel = \frac{N \varepsilon\, t_\text{device} V^2}{2 g}, \qquad
\bigl[ \mathbf{K}_e \bigr]_{\perp\perp} = -\frac{N \varepsilon\, t_\text{device}\, \ell\, V^2}{g^3}, \tag{112}
$$

along the fingers and across them. With the field method the computed
values exceed these by a few percent, the share of the fields at finger
tips and outer walls. The plates method returns exactly these when finger
tips are farther than $d_i$ from the wall they point at.

### Displaced bodies

Nothing above needs the bodies to be at rest. With every body at a
displacement from rest, stacked in $\mathbf{q}$, the same construction gives
a matrix $\mathbf{C}(\mathbf{q})$. A plate is evaluated by (55) at the
displacements of its two bodies. Sidewall outlines of body $B$ are carried
by the in-plane part (109) of its displacement,

$$
z \mapsto c_B + a + e^{\mathrm{i} a_\psi} (z - c_B), \qquad
a = a_x + \mathrm{i}\, a_y, \qquad
\begin{bmatrix} a_x \\ a_y \\ a_\psi \end{bmatrix} = \mathbf{J}_B\, \mathbf{q}_B, \tag{113}
$$

which is the motion (83) at finite size: travel, and rotation about the
center of mass, which travels with the body. Each cluster is then solved on
the carried outlines. Its gradient and second derivative there, taken about
the carried centers $c_B + a$, are the derivatives of its matrix in
$(a_x, a_y, a_\psi)$ at that displacement, because (83) composes: a further
small motion of a carried body is again (83) about its carried center. With
$\mathbf{J}_B$ constant, (110) holds at every $\mathbf{q}$, and so do (57)
and (58) for plates. Forces are therefore conservative,

$$
f_k(\mathbf{q}) = \frac{\partial}{\partial q_k}\, \frac{1}{2}\, \mathbf{v}^\top \mathbf{C}(\mathbf{q})\, \mathbf{v}, \qquad
\bigl[ \mathbf{K}_e(\mathbf{q}) \bigr]_{kl} = -\frac{\partial f_k}{\partial q_l}, \tag{114}
$$

exactly, for the function $\mathbf{C}(\mathbf{q})$ that the model computes.

What was found at rest stays. Clusters are those of (105) at rest. With the
field method the mesh is carried with the bodies and not made anew, and a
cluster solved in windows keeps the elements of each window, so
$\mathbf{C}(\mathbf{q})$ is smooth and (95) remains its exact derivative.
The plates method finds its pieces anew on the carried outlines, since they
are exact for any outlines; its kinks are those of (97). The limits are
those of the parts: a plate is exact in rise, of second order in tilt, and
of first order in in-plane travel, and the lever of (109) is that of a small
tilt. Electrodes of different bodies that touch after the motion are an
error.

!!! note "Implementation choices"
    The outline of a device electrode is that of the union of its body
    polygon and its stubs after a closing by the margin $\epsilon$ of (4):
    expansion by $\epsilon$ and erosion by the same, with mitered corners. A
    stub can meet its polygon in a single point, or along a line that
    rounding leaves a hair open, and the union alone would then keep a
    crack whose two walls nearly coincide. Plates take the union as it is,
    since a crack has no area. Regions of different electrodes never touch:
    a beam or a bridge lies between them, and its silicon belongs to no
    electrode. Where a beam met an electrode, the outline is closed by a
    short straight wall that carries charge like any other.

    The map (109) is linear. Terms of second order in tilt are left out of
    (110), as they are for plates, and so is any dependence of sidewall
    capacitance on travel out of the plane: facing sidewalls overlap fully
    at rest, so that dependence is even and exerts no force there.

    The second derivative is kept as one block per plate and per cluster,
    over the electrodes and degrees of freedom each involves, and (111) sums
    the blocks. A dense array over all pairs of electrodes and all pairs of
    degrees of freedom would not fit for a large structure.

    Left out of $\mathbf{C}$: sidewalls of the handle layer, sidewalls and
    faces of beams and bridges, coupling between clusters, and the field
    that leaves the plane of a layer at the top and bottom of a sidewall.

    With the field method a cluster of at most `max_elements` is one dense
    solution and a larger one is solved in windows. No cluster is left out:
    a matrix without it would look complete and not be. The number of
    windows and the asymmetry of (94) are kept for every cluster.

    The method is one choice for the whole structure and has no default. No
    rule in the library could tell a comb, which plates describe well, from
    a toothed motor, which they do not.

    Plates and sidewalls are independent of each other and of any
    mechanical solution: each can be left out, and each needs only its own
    parameters.

## Static equilibrium

With potentials on its nets a structure comes to rest where the pull of its
beams balances the electrostatic forces. Nothing flows at rest, so no link
carries a voltage and the electrodes of a net share its potential:
$\mathbf{v} = \mathbf{S}\, \mathbf{v}^\text{net}$ with $\mathbf{S}$ of (108).
Beams are linear, with the stiffness $\mathbf{K}^\text{full}$ of (35), and
the electrostatic force is that of (114) at the displaced bodies. With
mechanical loads $\mathbf{f}_0$ on the bodies, the displacement $\mathbf{q}$
solves

$$
\mathbf{K}^\text{full} \mathbf{q} = \mathbf{f}(\mathbf{q}) + \mathbf{f}_0, \qquad
f_k(\mathbf{q}) = \frac{1}{2}\, \mathbf{v}^\top \frac{\partial \mathbf{C}}{\partial q_k}(\mathbf{q})\, \mathbf{v}. \tag{115}
$$

The force grows as gaps close, so (115) is not linear and is solved by
Newton iteration from rest, $\mathbf{q}_0 = \mathbf{0}$:

$$
\bigl( \mathbf{K}^\text{full} + \mathbf{K}_e(\mathbf{q}_i) \bigr)\, \Delta \mathbf{q}
= \mathbf{f}(\mathbf{q}_i) + \mathbf{f}_0 - \mathbf{K}^\text{full} \mathbf{q}_i, \qquad
\mathbf{q}_{i+1} = \mathbf{q}_i + \Delta \mathbf{q}, \tag{116}
$$

with $\mathbf{K}_e$ of (114). Every iterate evaluates the capacitance model
at $\mathbf{q}_i$ for its gradient and second derivative, which are exact
there, so the iteration converges quadratically. Displacements mix meters
and radians, and steps are therefore measured by the strain energy they
would store in the beams. Iteration stops when

$$
\Delta \mathbf{q}^\top \mathbf{K}^\text{full} \Delta \mathbf{q}
\le \tau_s^2\; \mathbf{q}_{i+1}^\top \mathbf{K}^\text{full} \mathbf{q}_{i+1}, \tag{117}
$$

with $\tau_s$ the tolerance of the solve, or one step after the left side
has fallen below $\tau_s$ times the right side without the square. A step
that small is followed by one of about its square, which meets (117), unless
rounding in the forces keeps it larger, and then further steps would not
shrink it either.

An equilibrium holds only if a small push meets a restoring force in every
direction. The tangent stiffness $\mathbf{K}^\text{full} +
\mathbf{K}_e(\mathbf{q})$ answers that, and against the mass matrix of (34)
its smallest eigenvalue,

$$
\bigl( \mathbf{K}^\text{full} + \mathbf{K}_e(\mathbf{q}) \bigr)\, \boldsymbol{\phi}
= \lambda\, \mathbf{M}^\text{full} \boldsymbol{\phi}, \qquad
\lambda_\text{min} > 0, \tag{118}
$$

is the squared angular frequency of the softest motion under these
potentials, as in (38). The equilibrium is stable while it is positive.
Against the beam stiffness instead,

$$
\bigl( \mathbf{K}^\text{full} + \mathbf{K}_e(\mathbf{q}) \bigr)\, \boldsymbol{\psi}
= \mu\, \mathbf{K}^\text{full} \boldsymbol{\psi}, \qquad
\mu_\text{min} = \min_{\boldsymbol{\psi}} \frac{\boldsymbol{\psi}^\top \bigl( \mathbf{K}^\text{full} + \mathbf{K}_e \bigr) \boldsymbol{\psi}}{\boldsymbol{\psi}^\top \mathbf{K}^\text{full} \boldsymbol{\psi}}, \tag{119}
$$

the smallest eigenvalue is the share of beam stiffness left in the
direction the field softens most: 1 without potentials, 0 at pull-in, and
of the same sign as $\lambda_\text{min}$. The two answer different
questions. The motion that the field softens most can be stiff to begin
with, such as a plate moving toward the handle, and its frequency stays
far above that of the softest motion until close to pull-in, while
$\mu_\text{min}$ falls steadily. A square plate on four beams over the
handle at 95 percent of its pull-in voltage has $\mu_\text{min} = 0.44$,
while its softest motion, a travel in the plane that the field does not
touch, keeps the frequency it has without voltage.
For a plate of area $A$ on a spring $k$ over a gap $g_0$, (115) reads
$k\, z = -\varepsilon A V^2 / \bigl( 2 (g_0 + z)^2 \bigr)$ and (118) is
$k > \varepsilon A V^2 / (g_0 + z)^3$.

### Pull-in

Potentials are raised from zero in equal steps, $\mathbf{v}(s) = s\,
\mathbf{v}$ with the scale $s$ going from 0 to 1, and the equilibrium of one
step starts the iteration of the next, so the solution that is followed is
the one that grows out of rest. Pull-in is where it ends:

$$
s^* = \sup \bigl\{ s : \mathbf{q}(s) \text{ solves (115) at } \mathbf{v}(s),\ \text{continues from } \mathbf{q}(0),\ \lambda_\text{min}\bigl(\mathbf{q}(s)\bigr) > 0 \bigr\}. \tag{120}
$$

It ends in one of two ways. At a fold $\lambda_\text{min}$ falls to zero and
no equilibrium exists beyond: the plate above meets its fold at one third
of the gap. At a loss of symmetry the equilibrium goes on but no longer
holds against a push: a comb between two stators at one potential stays
centered while the spring across its fingers is used up by (112). The two
closed forms are

$$
V^* = \sqrt{\frac{8\, k\, g_0^3}{27\, \varepsilon A}} \ \text{ at } z = -\frac{g_0}{3}, \qquad
V^*_\perp = \sqrt{\frac{k_\perp\, g^3}{N \varepsilon\, t_\text{device}\, \ell}}. \tag{121}
$$

A step that finds no equilibrium, or finds one that is not stable, is
halved: the scale is bisected between the largest with a stable
equilibrium and the smallest without, each attempt starting from the
former, until the two are closer than the resolution asked for. The result
is the last stable equilibrium, its scale, and the scale just above it at
which stability was lost.

### Contact

Nothing in the model describes two bodies that touch, so an equilibrium is
accepted only while every gap between two bodies stays open. The clearance
of a displacement is the smallest ratio of such a gap to its width at rest,

$$
\chi(\mathbf{q}) = \min \Bigl\{
\min_{\text{plates}} \frac{g_0 + h_0 - |h_x|\, e_x - |h_y|\, e_y}{g_0},\;
\min_{(a, b)} \frac{\operatorname{dist}\bigl( E_a(\mathbf{q}), E_b(\mathbf{q}) \bigr)}{\operatorname{dist}(E_a, E_b)}
\Bigr\}. \tag{122}
$$

The first runs over plates between two bodies, with the rise $h_0$ and
slopes $h_x$, $h_y$ of (45) and the half widths $e_x$, $e_y$ of the bounding
box of the overlap about its reference point: its numerator is the gap (45)
at the corner of that box that comes closest. The second runs over pairs of
device electrodes of one cluster that move with different bodies and lie
within $d_i$ of each other at rest, with the outlines carried by (113), and
is the exact distance between their silicon. An equilibrium with
$\chi \le \chi_c$, the contact fraction, ends the rise of the potentials as
pull-in does, and the scale at which it happens is found by the same
bisection.

### Floating nets

A net that nothing drives keeps the charge it has. Split the nets into
driven ones, $d$, whose potentials are given, and floating ones, $f$, whose
charges $\mathbf{Q}_f$ are given. With the matrix between nets
$\mathbf{C}^\text{net}(\mathbf{q})$ of (108) in blocks, the charge on the
floating nets fixes their potentials at every displacement:

$$
\mathbf{C}^\text{net}_{ff}\, \mathbf{v}_f + \mathbf{C}^\text{net}_{fd}\, \mathbf{v}_d = \mathbf{Q}_f
\quad \Longrightarrow \quad
\mathbf{v}_f(\mathbf{q}) = \mathbf{C}_{ff}^{\text{net}\,-1} \bigl( \mathbf{Q}_f - \mathbf{C}^\text{net}_{fd}\, \mathbf{v}_d \bigr). \tag{123}
$$

The force keeps the form (115), now with these potentials. Its change with
displacement gains a part, because moving a body changes
$\mathbf{v}_f$: with $\boldsymbol{\Gamma}_f$ the rows for floating nets of
$\Gamma_{\nu k} = \sum_\mu \bigl( \partial C^\text{net}_{\nu\mu} / \partial q_k \bigr) v_\mu$,

$$
\mathbf{K}_t = \mathbf{K}^\text{full} + \mathbf{K}_e
+ \boldsymbol{\Gamma}_f^\top\, \mathbf{C}_{ff}^{\text{net}\,-1}\, \boldsymbol{\Gamma}_f \tag{124}
$$

replaces the tangent stiffness in (116), (118), and (119). The added part is never
negative: a body that moves toward a floating electrode lowers its
potential, which weakens the pull. A plate that holds a fixed charge and
has no other capacitance loses no stiffness at all, and never pulls in.
$\mathbf{C}^\text{net}_{ff}$ must be positive definite, which it is when
every floating net has capacitance to a driven one. Charges rise from zero
together with the potentials.

!!! note "Implementation choices"
    A Newton step that would close a gap between a plate and the handle,
    or bring electrodes of two bodies into contact, is halved until it does
    not, up to twelve times. A halved step never ends the iteration.

    An iteration that has not ended after `max_iterations` counts as no
    equilibrium at that scale, and so does one whose steps were halved or
    did not shrink four times in a row, which saves most of the cost of an
    attempt beyond pull-in. Beyond a fold that is the truth. It also
    happens with the plates method on teeth that stand tip to tip: the
    force jumps at the cusp of (97) and the iteration steps back and forth
    across it at any potential, so the result stops at a scale near zero
    without that being pull-in. Such a device needs the field method.

    The bounding box of a plate reaches at least as far as the plate, so
    (122) reports contact of a tilted plate no later than it happens, and
    earlier for a plate that does not fill its box. A tilted plate is also
    where the capacitance of (55), of second order in tilt, is least
    accurate: at a clearance of one half its integrand is 12 percent low
    at the closest corner.

    When bisection has closed in from below, the scale that failed is tried
    once more from the nearest stable equilibrium. An iteration can fail
    from a distant start where it succeeds from a near one, and a failure
    is only taken for pull-in when it comes from within the resolution.

    Forces of the field method carry rounding of its dense solutions, which
    grows with the size of a cluster. A structure with a very soft motion
    turns that into steps above a tight $\tau_s$, which is what the second way
    of stopping is for.

    Every iterate costs one evaluation of the capacitance model at
    displaced bodies. That is a fraction of a second with plates alone or
    with the plates method, and as long as the model at rest took with the
    field method.

    Beams stay linear at any displacement, and the limits of the
    capacitance model at displaced bodies carry over: tilt of a plate to
    second order, and in-plane travel of a plate to first.

    A body that no chain of beams joins to the anchor has no stiffness, and
    nothing balances a force on it. Such a body, for instance an island
    left on a neck too thin to qualify as a beam, is kept at rest and
    listed with the result. (115) to (119) are taken over the other bodies.

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
| $d_i$ | distance within which device electrodes are taken together for sidewalls; with the plates method also the widest gap between walls that count, with the field method the widest overlap of a window | `ElectrostaticParams.interaction_distance` |
| $n_w$ | most elements in the core of a window | `ElectrostaticParams.window_elements` |
| $\tau$ | change of the charge on a core below which the overlap of its window stops growing | `ElectrostaticParams.window_tolerance` |
| $\tau_s$ | relative size of the last Newton step at which static equilibrium is accepted | `StaticOptions.tolerance` |
| $\chi_c$ | fraction of its width at rest to which a gap may close before an equilibrium counts as contact | `StaticOptions.contact_fraction` |
| $E$, $G$ | Young's and shear modulus of the device layer | `Material.youngs_modulus`, `Material.shear_modulus` or `Material.poisson_ratio` |
| $J$ | torsion constant model | `Material.torsion_model` |
| $\beta$ | short-beam torsional stiffening coefficient | `Material.torsion_beta` |
