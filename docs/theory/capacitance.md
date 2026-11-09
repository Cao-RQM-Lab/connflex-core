# Capacitance of a structure

This page assembles the capacitance matrix of a structure from its plates,
its sidewalls found cluster by cluster, and its fringe, gives the forces and
the electrostatic stiffness that its derivatives in the motion of the bodies
produce, and defines its evaluation with the bodies displaced from rest, on
which every coupled analysis rests. `Structure.capacitance` builds the model
at rest and `CapacitanceModel.at` evaluates it elsewhere.

The electrodes of a structure, $n$ in all, are the nodes of one capacitance
matrix $\mathbf{C}$, defined by the charges $\mathbf{Q} = \mathbf{C}
\mathbf{v}$ that potentials $\mathbf{v}$ put on them. It is the sum of the
plates between device and handle electrodes and of the sidewalls between
device electrodes.

A plate of capacitance $c$ between electrodes $a$ and $b$ stores energy
$\tfrac{1}{2} c\, (v_a - v_b)^2$, so with $\mathbf{e}_a$ the unit vector of
electrode $a$ it adds

$$
c\, (\mathbf{e}_a - \mathbf{e}_b)(\mathbf{e}_a - \mathbf{e}_b)^\top \tag{109}
$$

to $\mathbf{C}$: $c$ to the two diagonal entries and $-c$ to the two
others. Each plate enters twice, with $C_\text{ox}$ of (43) into
$\mathbf{C}^\text{ox}$ and with $\varepsilon A / g_0$ into
$\mathbf{C}^\text{rel}$.

Sidewalls are found in the plane of the device layer. In that plane the
field couples every conductor to every other, however far, which the handle beneath and the space above prevent in
the real stack. Electrodes are therefore taken together only when their
silicon is close. Let $d_i$ be the interaction distance and join device
electrodes $a$ and $b$ when

$$
\operatorname{dist}(E_a, E_b) \le d_i. \tag{110}
$$

A cluster is a connected set of this relation with two or more electrodes:
electrodes joined directly or through a chain of neighbors. Each cluster is
one plane problem with its electrodes as conductors and their nets as the
nets of (105), and its matrix $\mathbf{c}$ is that of (108), each pair of
nets found by its own method. The same $d_i$ is the widest gap across which
two walls count as plates and the widest overlap of a window. The matrix
times the thickness of the device layer fills the rows and columns of its
electrodes,

$$
\bigl[ \mathbf{C}^\text{sw} \bigr]_{ab} = t_\text{device}\, c_{ab} \quad \text{for } a, b \text{ in one cluster}, \tag{111}
$$

and electrodes of different clusters share no sidewall capacitance. With
$\mathbf{C}^\text{fr}$ the fringe above and below the device layer of
(163), found for each cluster and for each device electrode alone over the
handle, zero when it is left out,

$$
\mathbf{C} = \mathbf{C}^\text{ox} + \mathbf{C}^\text{rel} + \mathbf{C}^\text{sw} + \mathbf{C}^\text{fr}, \qquad
\mathbf{C}^\top = \mathbf{C}, \qquad \mathbf{C}\, \mathbf{1} = \mathbf{0}. \tag{112}
$$

Rows sum to zero because (109), (71), (98), and (163) do: a potential common to all
electrodes moves no charge, nothing outside the structure holds charge, and
$-C_{ab}$ is the capacitance between electrodes $a$ and $b$.

Electrodes of one net are joined by links and share a potential at rest.
With $S_{a\nu} = 1$ when electrode $a$ belongs to net $\nu$ and zero
otherwise, the capacitance matrix between nets is

$$
\mathbf{C}^\text{net} = \mathbf{S}^\top \mathbf{C}\, \mathbf{S}, \tag{113}
$$

which adds the entries of $\mathbf{C}$ over the electrodes of each pair of
nets. Capacitance between two electrodes of one net drops out.

### Forces and stiffness

$\mathbf{C}$ depends on the positions of the bodies. Let $\mathbf{q}$ stack
the degrees of freedom (12) of every movable body, in the order of the
assembled matrices (34) and (35). A plate between two bodies contributes
its gradient and Hessian (58) in the pattern of (109). A cluster contributes
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
\end{bmatrix}, \qquad d_z = z_m - c_{B,z}, \tag{114}
$$

so that with the gradient and second derivative of the cluster, (81) and
(84) for the field method or those of (103) for plates, taken about the
centers of mass,

$$
\frac{\partial \mathbf{C}^\text{sw}}{\partial \mathbf{q}_B} = t_\text{device}\, \frac{\partial \mathbf{c}}{\partial \mathbf{a}_B}\, \mathbf{J}_B, \qquad
\frac{\partial^2 \mathbf{C}^\text{sw}}{\partial \mathbf{q}_B\, \partial \mathbf{q}_{B'}^\top}
= t_\text{device}\, \mathbf{J}_B^\top \frac{\partial^2 \mathbf{c}}{\partial \mathbf{a}_B\, \partial \mathbf{a}_{B'}^\top}\, \mathbf{J}_{B'}. \tag{115}
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
\bigl[ \mathbf{K}_e \bigr]_{kl} = -\frac{1}{2}\, \mathbf{v}^\top \frac{\partial^2 \mathbf{C}}{\partial q_k\, \partial q_l}\, \mathbf{v}. \tag{116}
$$

$\mathbf{K}_e$ adds to $\mathbf{K}^\text{full}$ of (35). It is symmetric,
and a negative diagonal entry means that the field softens that motion. For
a comb of $N$ gaps of width $g$ and overlap $\ell$ between a stator at $V$
and a shuttle at zero, the ideal values that (116) approaches are

$$
f_\parallel = \frac{N \varepsilon\, t_\text{device} V^2}{2 g}, \qquad
\bigl[ \mathbf{K}_e \bigr]_{\perp\perp} = -\frac{N \varepsilon\, t_\text{device}\, \ell\, V^2}{g^3}, \tag{117}
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
by the in-plane part (114) of its displacement,

$$
z \mapsto c_B + a + e^{\mathrm{i} a_\psi} (z - c_B), \qquad
a = a_x + \mathrm{i}\, a_y, \qquad
\begin{bmatrix} a_x \\ a_y \\ a_\psi \end{bmatrix} = \mathbf{J}_B\, \mathbf{q}_B, \tag{118}
$$

which is the motion (83) at finite size: travel, and rotation about the
center of mass, which travels with the body. Each cluster is then solved on
the carried outlines. Its gradient and second derivative there, taken about
the carried centers $c_B + a$, are the derivatives of its matrix in
$(a_x, a_y, a_\psi)$ at that displacement, because (83) composes: a further
small motion of a carried body is again (83) about its carried center. With
$\mathbf{J}_B$ constant, (115) holds at every $\mathbf{q}$, and so do (57)
and (58) for plates. Forces are therefore conservative,

$$
f_k(\mathbf{q}) = \frac{\partial}{\partial q_k}\, \frac{1}{2}\, \mathbf{v}^\top \mathbf{C}(\mathbf{q})\, \mathbf{v}, \qquad
\bigl[ \mathbf{K}_e(\mathbf{q}) \bigr]_{kl} = -\frac{\partial f_k}{\partial q_l}, \tag{119}
$$

exactly, for the function $\mathbf{C}(\mathbf{q})$ that the model computes.

What was found at rest stays. Clusters are those of (110) at rest. With the
field method the mesh is carried with the bodies and not made anew, and a
cluster solved in windows keeps its cores, each window its elements while
(164) holds, and each row its cells while (165) holds, so
$\mathbf{C}(\mathbf{q})$ is smooth between changes of set and (95) its
exact derivative there. With `include_fringe` it also steps where a gap
crosses a line of the grid of (160), by the residual of the rounding, of
relative order $(\Delta / g)^2$ of the fringe part; Newton iteration can
chatter at such a line, and the energy ledger of a transient counts the
steps.
The plates method finds its pieces anew on the carried outlines, since they
are exact for any outlines; its kinks are those of (97). The limits are
those of the parts: a plate is exact in rise, of second order in tilt, and
of first order in in-plane travel, and the lever of (114) is that of a small
tilt. Electrodes of different bodies in one cluster that touch after the
motion, and a plate whose gap closes at a corner of its overlap, are
errors; electrodes farther apart than $d_i$ at rest are in no cluster and
are not checked.

!!! note "Outline of an electrode"
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

!!! note "Lever of the mid-plane and tilt"
    The map (114) is linear. Terms of second order in tilt are left out of
    (115), as they are for plates. Sidewall capacitance in the layout plane
    does not depend on travel out of it: rise and tilt enter through the
    fringe of (157) to (163), which gives the lift of fingers over the
    handle and the vertical spring of a comb.

!!! note "Hessian kept in blocks"
    The second derivative is kept as one block per plate and per cluster,
    over the electrodes and degrees of freedom each involves, and (116) sums
    the blocks. A dense array over all pairs of electrodes and all pairs of
    degrees of freedom would not fit for a large structure.

!!! note "What the matrix leaves out"
    Left out of $\mathbf{C}$: sidewalls of the handle layer, sidewalls and
    faces of beams and bridges, and coupling between clusters; without
    `include_fringe`, also the field that leaves the plane of the device
    layer at the top and bottom of a sidewall.

!!! note "No cluster is left out"
    No cluster is left out, whatever its size: a matrix without it would
    look complete and not be. The feature size and method of every pair of
    nets, the number of windows, and the asymmetry of (94) are kept for
    every cluster.

!!! note "Plates and sidewalls are independent"
    Plates and sidewalls are independent of each other and of any
    mechanical solution: each can be left out, and each needs only its own
    parameters.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $\mathbf{C}$ | capacitance matrix between electrodes | `CapacitanceModel.matrix` |
| $\mathbf{C}^\text{ox}$, $\mathbf{C}^\text{rel}$, $\mathbf{C}^\text{sw}$, $\mathbf{C}^\text{fr}$ | its parts: oxide, released gap, sidewalls, fringe | `CapacitanceModel.oxide`, `released`, `sidewall`, `fringe` |
| $d_i$ | distance within which device electrodes form a cluster | `ElectrostaticParams.interaction_distance` |
| $\mathbf{C}^\text{net}$, $\mathbf{S}$ | matrix between nets and the net of each electrode | `CapacitanceModel.net_matrix`, `CapacitanceModel.net` |
| $\mathbf{q}$ | degrees of freedom of every movable body | `CapacitanceModel.dofs`, the argument of `CapacitanceModel.at` |
| $\mathbf{J}_B$, $d_z$ | lever from a body's degrees of freedom to its motion at the device mid-plane | internal to `electrostatics.capacitance` |
| $\partial \mathbf{C} / \partial \mathbf{q}$, $\partial^2 \mathbf{C} / \partial \mathbf{q}^2$ | gradient and Hessian of the matrix | `CapacitanceModel.gradient`, `CapacitanceModel.blocks` |
| $f_k$, $\mathbf{K}_e$ | generalized force and electrostatic stiffness | `CapacitanceModel.force`, `CapacitanceModel.stiffness`, `CapacitanceModel.charge_gradient` |
| $\mathbf{C}(\mathbf{q})$ | the model at displaced bodies | `CapacitanceModel.at`, `CapacitanceModel.displacement` |
