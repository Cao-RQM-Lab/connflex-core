# Plates between device and handle

Every piece of device silicon that lies over handle silicon is a capacitor
with it, through buried oxide where oxide survives and across the released gap
elsewhere. This page gives the capacitance of such a plate and its derivatives
in the six degrees of freedom of the two bodies it joins, exact in vertical
travel, to second order in tilt, and to first order in travel in the plane.
`Structure.plates` finds the plates of a structure and `Plate.evaluate` returns
capacitance, gradient, and Hessian at a displacement of its two bodies; the
[capacitance model](capacitance.md) sums them.

Silicon is treated as a perfect conductor, so an electrode
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

!!! note "Moments without building the overlap"
    The moments of $G_{ab}$ are found without building it. The sums of (11b)
    are evaluated for $U_{ab}$, for the oxide polygons inside it, and for the
    release openings inside it, and the last two are subtracted. Oxide and
    openings are disjoint, since (2) keeps oxide at least $d$ away from every
    opening, and each set is disjoint in itself. A polygon that $U_{ab}$
    covers is used as drawn and one that crosses its boundary is clipped to
    it first. A plate with a hundred thousand openings therefore costs one
    pass over their vertices, as for mass.

!!! note "Building the handle outline under the plate"
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

!!! note "Region of a device electrode"
    The region of a device electrode is the union of its body polygon and
    its stubs. Silicon of beams and bridges belongs to no electrode, so its
    capacitance to the handle is not counted.

!!! note "Plates within one body"
    When both electrodes move with one body, (46) gives $\mathbf{s} = 0$ for
    every motion and the plate is a fixed capacitor. A bonded pad is the
    common case: oxide under its middle and the undercut rim of width $d$
    around it, both between silicon of the anchor.

!!! note "No fringing field in the plate"
    The plate itself holds no fringing field. With `include_fringe` the
    field at device walls over the handle is added by the open cuts of
    (162); the field into release openings is not modeled, and openings not
    much wider than the gap recover part of the area they remove, so the
    released part of (43) is a lower bound there.

!!! note "Reference point and mid-height"
    The reference point is the center of the bounding box of the overlap
    rather than the origin or a center of mass, so that $\mathbf{S}$,
    $\mathbf{I}$, and the edge integrals stay of the size of the plate and
    the products in (58) do not cancel. In-plane travel is taken at
    mid-height of the gap for both bodies: taking each at its own face would
    turn a rigid tilt common to the two into travel of $t_\text{box}$ times
    that tilt. A rise that closes the gap at the reference point,
    $\bar{g} \le 0$, is an error.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $E_a$, $E_b$ | device and handle electrode of a plate | `Plate.electrodes` |
| $U_{ab}$, $O_{ab}$, $G_{ab}$ | overlap and its bonded and released parts | `Plate.oxide` and `Plate.released` hold their capacitances |
| $g_0$ | thickness of the buried oxide, the gap at rest | `Layer.thickness` of the box layer |
| $\varepsilon_\text{ox}$, $\varepsilon_a$ | relative permittivity of the oxide and of the medium in the gap | `Material.relative_permittivity` of the box layer, `ElectrostaticParams.ambient_permittivity` |
| $\mathbf{s} = (a_x, a_y, a_\psi, h_0, h_x, h_y)$ | local coordinates: travel and turn in the plane, rise, and slopes of the gap | internal to `Plate.evaluate`; rise and slopes enter the clearance of (126) |
| $\mathbf{J}(\mathbf{c})$ | map from the degrees of freedom of a body to the local coordinates | `mechanics.kinematics.point_jacobian` |
| $A$, $\mathbf{S}$, $\mathbf{I}$ | area moments of the released overlap about the reference point | `mechanics.inertia.AreaMoments` |
| $\ell_k$, $\boldsymbol{\lambda}_k$ | edge integrals over the handle outline under the plate | internal to `electrostatics.plates` |
| $V_{ab}$ | vertices of the convex hull of the overlap, where a tilted gap is narrowest | `Plate.corners` |
