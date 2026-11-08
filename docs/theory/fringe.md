# Fringe above and below the device layer

The layout-plane field treats every sidewall as if it went on without end in
height. Real walls end at the top and bottom faces of the device layer, and
the field that leaves those ends arches across gaps and reaches down to the
handle. This page adds that field: a two-dimensional cut across every wall,
solved in the plane of its normal and the vertical, what the cut holds beyond
the layout-plane field and the plates, and how the cuts are summed along the
walls, expanded from a grid, and rounded where walls pass each other.
`include_fringe` of `ElectrostaticParams` switches it on, and its part of the
matrix is `CapacitanceModel.fringe`.

The layout-plane field of (61) is the field of walls without end in height: it
puts $\varepsilon t / g$ per unit length between walls a gap $g$ apart, and
nothing above or below the device layer. Real walls end at its top and bottom
faces. Field leaves the top face of one wall, arches over the gap, and lands on
the top of the other; at the bottom it divides between the other wall and the
handle a box gap $g_0$ below. A thinner layer or a wider gap leaves more of the
field to these ends, and since the handle lies under the device layer alone,
it draws the field of the lower ends down and lifts what lies between fingers
at another potential.

Near a wall the field is taken as two-dimensional in the plane of its normal
and the vertical, a cut, and what the field of the cut holds beyond the
layout-plane field and the plates is added along the wall. At a point
$\mathbf{p}$ of a wall of conductor $w$ with outward normal $\hat{\mathbf{n}}$,
the line $\mathbf{p} + s\, \hat{\mathbf{n}}$ for $|s| \le d_i$ crosses silicon
in intervals $[s_k^-, s_k^+]$ of conductors $c_k$, the wall's own ending at
$s = 0$. In the plane of $s$ and height $z$ each interval is a block, and the
handle, where it lies beneath the cut, a slab:

$$
B_k = [s_k^-, s_k^+] \times [z_{c_k},\, z_{c_k} + t], \qquad
H = [-(d_i + 4 t),\, d_i + 4 t] \times [z_h - 2 g_0,\, z_h - g_0], \tag{157}
$$

with $z_c$ the rise of the body of conductor $c$ at the cut and $z_h$ that of
the handle: a closed slab one box gap tall, whose lower face and ends carry
charge too. The cut is solved as (63) on the outlines of these blocks and of
the slab, which reaches $d_i + 4 t$ past the cut on either side. Its core $K$ is
the wall face $s = 0$ and the halves of the top and bottom faces of the wall's
block nearest it, and the charge on the core per unit length per volt on
conductor $j$, as (92), is

$$
q_j = \mathbf{e}_K^\top \mathbf{M}\, \mathbf{e}_j. \tag{158}
$$

The layout-plane field already holds $t / d$ between the wall and $f$, the
first conductor ahead of it at $s = d$, and a plate (42) holds $(w/2) /
(g_0 + z_w - z_h)$ under the half face of the core, $w$ the width of the wall's
block. What the cut adds, per unit length and permittivity, is the rest:

$$
\epsilon_j = q_j - \frac{t}{d}\, (\delta_{jw} - \delta_{jf}) - \frac{w / 2}{g_0 + z_w - z_h}\, (\delta_{jw} - \delta_{jh}). \tag{159}
$$

Its derivatives in the travel of each body along $s$ and upward are those of
(76) and (84) on the outlines of the cut, and of the two subtractions. The
handle is a plane: travel along the cut does not move it, and the slab, which
is centered on the wall, goes with the wall.

Walls repeat, and so do cuts. A cut is rounded to a grid of $\Delta = \min(t,
g_0, g_\text{min}) / 4$, with $g_\text{min}$ the narrowest gap between
electrodes of different bodies of the group at rest, and solved once for
every rounded cut met; a cut that rounding would close, silicon or a gap
narrower than $\Delta$, keeps its own geometry. The expansion (160) carries
the rounding of the nearest end of each body to first order, so what the
grid costs the capacitance is of second order in $\Delta$; the force, which
also feels the far ends and the ends of other bodies that are rounded
without correction, pays between first and second order, and the stiffness,
whose second derivative is that of the rounded cut, first order. On a comb
of fingers 4 micrometers wide at gaps of 3 micrometers, 25 micrometers thick
over a box of 2 micrometers, the grid costs $2 \times 10^{-5}$ of the fringe
matrix, $2 \times 10^{-3}$ of the force, and $2 \times 10^{-3}$ of the
stiffness at displacements off the grid; a grid four times finer costs
$2 \times 10^{-6}$, $4 \times 10^{-4}$, and $5 \times 10^{-4}$. What the
grid saves is bounded, since cuts of one layout differ in how many spans
they cross and of which conductors far more than in how those spans round.
Each body $b$ of a cut has the end of its silicon nearest the wall a distance
$\delta_b$ from where rounding put it, measured from the wall, and rises by
$z_b$. With $\mathbf{x} = (\delta_b, z_b)$ over the bodies of the cut,

$$
\epsilon_j(\mathbf{x}) = \epsilon_j^0 + \mathbf{g}_j^\top \mathbf{x} + \tfrac{1}{2}\, \mathbf{x}^\top \mathbf{H}_j\, \mathbf{x}, \tag{160}
$$

with $\epsilon^0$, $\mathbf{g}$, and $\mathbf{H}$ those of the rounded cut. A
body rises at the point where the cut stands on its wall, carried back to
rest, by $z_b = u_z + \theta_x (y - y_b) - \theta_y (x - x_b)$ with $(x_b,
y_b)$ the center of body $b$, as a plate does at its reference point (46).

Every piece of wall that faces another conductor, (99) to (101), is cut at its
middle and counts with its length $L$; the offset of the body it faces is the
gap of (102) at that middle, measured to the end of its silicon ahead of the
wall even where that body has silicon nearer behind it, and those of other
bodies their travel along the cut, as if their walls were square to it. A piece whose ends are ends of its
own edge, from 0 to $a$ along it, and of the edge it faces, from $b_1$ to
$b_2$, has its length rounded over $\lambda = g_0$ where ends pass,

$$
L = r^+(a - b_1) - r^\circ(a - b_2) - r^\circ(-b_1) - b_2 + r^+(b_2), \tag{161}
$$

$$
r^+(x) = \begin{cases} 0 & x \le 0 \\ \frac{2 x^2}{\lambda} - \frac{x^3}{\lambda^2} & 0 < x < \lambda \\ x & x \ge \lambda \end{cases}, \qquad
r^\circ(x) = \begin{cases} 0 & x \le -\frac{\lambda}{2} \\ \frac{(x + \lambda / 2)^2}{2 \lambda} & |x| < \frac{\lambda}{2} \\ x & x \ge \frac{\lambda}{2} \end{cases},
$$

which is the overlap of the two away from where ends pass, has a continuous
slope everywhere, and meets zero with zero slope where the walls stop facing
each other: $r^+$ rounds where walls begin or stop facing, inside the overlap,
and $r^\circ$ where two ends on one side pass, evenly either way. A tooth that
slides past another so pulls without a jump.

A wall that faces nothing, or the part of an edge that no piece covers, is
cut behind its wall alone, with nothing ahead and no layout-plane part to
take out, wherever the handle lies beneath the middle of its edge. Every edge
$e$ of length $\ell_e$ counts its open cut $\epsilon^\circ_e$, and every piece
$k$ on it counts its own cut less the open one,

$$
\hat{c}_{wj} = \sum_{e} \ell_e\, \epsilon^\circ_{ej} + \sum_{k} L_k \bigl(\epsilon_{kj} - \epsilon^\circ_{e(k) j}\bigr), \tag{162}
$$

over the edges and pieces of conductor $w$, each cut expanded by (160). Like
(93) the counts are not symmetric, and the capacitance added is

$$
C_{ab} = \varepsilon\, \frac{\hat{c}_{ab} + \hat{c}_{ba}}{2}, \qquad
C_{ah} = C_{ha} = \varepsilon\, \hat{c}_{ah}, \qquad
C_{aa} = -\sum_{b \ne a} C_{ab}, \tag{163}
$$

between device electrodes $a$ and $b$ and to a handle electrode $h$. Its
gradient and Hessian follow (162) and (163) through the lengths of (161) and
the expansions of (160), in the six degrees of freedom of every body.

Two walls 4 micrometers apart, the wall's block 30 micrometers wide, as far as
the cut reaches with $d_i$ of 30 micrometers, and the block ahead 22, show how
the correction goes with thickness:

| $t$ | Excess to the wall ahead | Against $t / g$ | With the handle 2 micrometers below: to the wall ahead | To the handle |
| --- | --- | --- | --- | --- |
| 10 um | $1.85\, \varepsilon$ | 0.74 | $0.73\, \varepsilon$ | $0.98\, \varepsilon$ |
| 20 um | $1.85\, \varepsilon$ | 0.37 | $0.73\, \varepsilon$ | $0.96\, \varepsilon$ |
| 40 um | $1.85\, \varepsilon$ | 0.18 | $0.75\, \varepsilon$ | $0.92\, \varepsilon$ |
| 80 um | $1.85\, \varepsilon$ | 0.09 | $0.77\, \varepsilon$ | $0.88\, \varepsilon$ |

The fringe at the two ends of a wall does not grow with its height, so it
matters less as $t / g$ grows, and the handle takes well over half of it. On a
comb of fingers 4 micrometers wide with gaps of 3 micrometers, 25
micrometers thick, the fringe adds 14 percent to the pull along the fingers
without a handle. Over a grounded handle, with the shuttle grounded and the
stator driven, it takes 20 percent off the pull instead: walls of the stator
that no finger faces couple to the handle beside them, and an engaging finger
shields them from it.

!!! note "Mesh of a cross-section"
    The cross-section holds air: buried oxide under anchored electrodes is
    left to the plates, since the release undercut leaves air under the
    edges of walls, where the fringe lies, wherever the release distance
    exceeds the box. Elements are at most $t / 8$ long on walls and on faces
    away from features, at most the smaller of $g_0$ and $t / 8$ on faces
    where they pass the wall or the edge of a block, growing by at most $e$
    per element in between, and shorten toward corners by `corner_grading`.
    The slab is meshed as finely under every block edge, so its mesh follows
    the blocks of a cut; the derivatives of a cut hold it where it is, and
    (160) expands from cuts on the grid, where it stands still. Between 4 and
    32 elements over the height the excess changes by $2 \times 10^{-4}$
    without the handle and by $3 \times 10^{-3}$ with it. The slab reaches
    $d_i + 4 t$ either side of the wall: its top beyond the silicon the cut
    knows takes charge from the top face of the wall's block, which the cut
    truncates at $d_i$, so the excess settles slowly as the slab grows. From
    $4 t$ to $8 t$ the excess to the wall ahead falls by $6 \times 10^{-3}$
    and that to the handle rises by $1 \times 10^{-2}$, and to $16 t$ by
    $8 \times 10^{-3}$ and $1.3 \times 10^{-2}$; a wall's block wider than
    about $t$ so keeps a fringe to its faced conductor that grows a little
    with $d_i$.

!!! note "What a cut takes as uniform"
    A cut is uniform along its wall. Corners of outlines add no field of
    their own, and the rounding of (161) stands in for it where ends pass.
    Walls that face each other at an angle steeper than plates take,
    (101), count as open. An open wall with no handle beneath counts nothing:
    its fringe reaches only conductors behind it within $d_i$. With
    `include_plates` off the cuts have no handle. Device electrodes with no
    neighbor within $d_i$ count their edges over the handle too.

!!! note "Exactness of the derivatives"
    Gradients are the exact derivatives of (163). Second derivatives are
    exact in travel and rise for the wall's body and the body it faces;
    walls of other bodies in a cut move as if square to it, which leaves
    those in turn a few tenths of a percent of the fringe's off. A block in
    a cut rises as a whole: the slope of a tilted body within the plane of
    the cut is left out, a term of first order in the tilt $h$ and of
    relative size $h w / (2 g_0)$ of the fringe to the handle for a block
    of width $w$.

!!! note "Silicon along a cut"
    The silicon a cut line crosses is read off the edges of the polygons,
    held in a spatial index, rather than from an overlay of the line with
    each polygon, which costs far more on outlines with many vertices. An
    edge whose envelope meets the line's is tested for a crossing by the
    half-open rule, the line taken to lie a hair to its left so that a
    vertex on it counts once; the state of the line just before its start
    is read from the polygons a hair before and to the left of that point;
    an edge along the line is silicon, as the closed polygon has it; and
    envelopes are widened by the tolerance so that an edge a line ends on,
    to rounding, is found. The stretches are those of the overlay.

!!! note "Cost and reuse"
    A cross-section takes tens of milliseconds. Cuts that repeat, as in
    combs, are solved once for a structure, clusters sharing them. Displaced
    bodies take the cross-sections whose rounding they share and expand from
    them by (160); a surrogate, (149), keeps the fringe of each of its
    solutions and expands it with them.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $t$, $g_0$ | thickness of the device layer and of the buried oxide | `Layer.thickness` of the device and box layers |
| $d_i$ | reach of a cut either way from its wall | `ElectrostaticParams.interaction_distance` |
| $B_k$, $H$ | blocks of silicon in a cut and the handle slab | `electrostatics.section`; the slab reaches `HANDLE_MARGIN` thicknesses past the cut |
| $q_j$, $\epsilon_j$ | charge on the core of a cut and its excess over plane field and plate | `Section` results, `Cluster.sections` counts distinct cuts |
| $\Delta$ | grid the cuts are rounded to | `electrostatics.fringe.KEY_STEPS` divides the smallest of $t$, $g_0$, and the narrowest gap |
| $\mathbf{x} = (\delta_b, z_b)$ | offset of each body from the rounded cut and its rise | internal to `electrostatics.fringe` |
| $L$, $r^+$, $r^\circ$ | rounded length of a piece and the rounding functions | internal to `electrostatics.fringe` |
| $\hat{c}_{wj}$, $C_{ab}$ | fringe counted along the walls and the capacitance added | `CapacitanceModel.fringe` |
