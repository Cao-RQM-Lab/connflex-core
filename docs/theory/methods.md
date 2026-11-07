# Facing plates and the method for each pair of nets

A field solution is the most accurate way to find sidewall capacitance and the
most expensive. In a comb, though, most of the capacitance sits in narrow gaps
between long walls, where the field is that of a parallel plate and a field
solution spends its effort on what a closed form already gives. This page
describes the plates method, which takes every wall that faces another
conductor as a plate with it, and the feature size by which each pair of
interacting nets is sent to one dense field solution, to windows, or to
plates. `sidewall_capacitance` and `Structure.capacitance` both choose this
way, and each pair's feature size and method are reported with the result.

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
of tooth tops. Such a device needs the field method, which the feature
size of (105) selects for it: tooth tops a few gaps wide are far below any
plates cutoff that keeps plates accurate.

Teeth that stand tip to tip at rest are the worst case. By plates their
capacitance is a cusp in travel: it falls to either side of rest, the two
one-sided forces are large and opposite, and the force at rest is whichever
of them the rounding of coordinates selects, and it can exceed the force
of the field many times over. The electrostatic spring there is the mean of
the two sides and is far less affected.

!!! note "Finding the edges in view"
    Edges in view of a wall are found in a spatial index of all edges,
    queried with the rectangle that the wall sweeps along its normal over
    $d_i$. An edge that lies along the normal, with end points less than
    $10^{-9} L_i$ apart in $u$, has no height (100) and is skipped. Marks
    closer than $10^{-9} L_i$ count as one.

!!! note "Series of the piece integral"
    Below $|x| = 1/4$ the second factor of (102) is summed from its series,
    so nothing cancels as $m_j \to 0$.

!!! note "Walls along a normal and coincident marks"
    Where a wall lies exactly along the normal of another, as the end face
    of a finger does beside the wall its flank faces, the slope limit of
    (101) keeps the end face out of view until a turn brings its slope seen
    from the other wall below $1/10$, so (97) is smooth in rotation at the
    square position and changes by a jump where that slope crosses $1/10$,
    of up to $\tfrac{1}{2} \varepsilon\, \theta \ln (1 + w / \gamma)$ for an
    end face of width $w$ beside a gap $\gamma$ at the turn $\theta$, less
    where $d_i$ cuts the view, small against the $\varepsilon L r /
    \gamma^2$ of a flank of length $L$ at distance $r$ from the pivot. Where
    marks set by end points of different bodies
    coincide, as when a finger tip lines up exactly with a corner across
    the gap, there is a kink in travel too. Of coincident marks the one that
    sets an end of a piece is an end point of the edge faced if there is
    one, else an end of the wall itself, else an end point of the edge
    faced on the other side of the mark.

!!! note "Why the slope limit"
    The limit on slope makes a wall that is square to another only up to
    the rounding of its coordinates behave as one that is square exactly.
    Taken as a plate it would be a sliver a few picometers wide, whose
    second derivative grows with its slope without bound, so on teeth turned
    slightly out of the grid such slivers would give an electrostatic spring
    of either sign and of any size.

!!! note "The cutoff is sharp"
    The cutoff is sharp. A wall parallel to another just inside $d_i$
    counts in full and just outside not at all, so $d_i$ belongs well above
    every gap that matters and below the distance to silicon that should
    not count.

!!! note "No mesh, pieces kept"
    The method finds no charge density and needs no mesh. The pieces of
    wall, with the conductor each faces and the gap at both ends, are kept
    for review.

## Method for each pair of nets

The three methods differ most where fringe field matters. Plates miss it
entirely; windows miss only the far field beyond their overlap; one dense
system misses nothing but grows with the cube of its size. Fringe field
matters least between long walls across narrow gaps, which is also where
outlines are longest and a field solution costs most. Each pair of nets
that interact therefore takes its method from the shape of its facing walls.

### Feature size

The pieces of (99) to (101) that continue each other form runs. A piece
continues another when it starts where the other ends, on the outline of
the same conductor, facing the same net, and when the points at which the
normals of the two meet the faced wall at their junction lie within half a
gap of each other. A run follows a wall around a drawn arc or a corner of
its own outline as long as the wall it faces runs on, and ends where the
gap jumps, as at the edge of a tooth. Let run $r$ have length $L_r$ and

$$
w_r = \int_r \frac{ds}{\gamma}, \tag{104}
$$

its capacitance per permittivity as a plate, which is also its length over
its mean gap, $L_r / \bar{\gamma}_r$ with $1 / \bar{\gamma}_r$ the mean of
$1 / \gamma$ along it. The feature size of nets $A$ and $B$ is the mean of
$w_r$ over the runs of $A$ facing $B$ and of $B$ facing $A$, weighted by
$w_r$:

$$
\lambda_{AB} = \frac{\sum_r w_r^2}{\sum_r w_r}. \tag{105}
$$

It is the length of their facing walls in units of the gap between them,
counted where the capacitance is: for a comb the finger overlap over the
finger gap, for a plate actuator its length over its gap, for teeth the
tooth width over the gap at the tips. Two nets whose silicon lies within
$d_i$ but whose walls do not face have $\lambda_{AB} = 0$. Electrodes of one
net are a pair $A = B$ like any other.

### Choice and assembly

With the dense cutoff $\lambda_d$ and the plates cutoff $\lambda_p$, the
pair takes

$$
m_{AB} =
\begin{cases}
\text{one dense system} & \lambda_{AB} < \lambda_d, \\
\text{windows} & \lambda_d \le \lambda_{AB} < \lambda_p, \\
\text{plates} & \lambda_p \le \lambda_{AB},
\end{cases} \tag{106}
$$

and every pair of conductors $a$, $b$ of a cluster, one in net $A$ and one in
net $B$, within $d_i$ of each other, takes the method of $(A, B)$. Plates
pairs take $c_{ab}$ of (97) and its derivatives (103). Pairs solved by the
field count charge on the elements of their two conductors. The elements of
the conductors of all dense pairs form the core of one window, (91) with a
single core, and those of all windowed pairs are cut into cores of at most
$n_w$ elements. A window is solved only when it holds both conductors of
one of its pairs and its core lies on one of them, so where only pairs taken
as plates meet, nothing is solved. Each gives $\hat{\mathbf{c}}$ of (93)
over all conductors. A pair takes

$$
p_{ab} = -\tfrac{1}{2} \bigl( \hat{c}_{ab} + \hat{c}_{ba} \bigr) \tag{107}
$$

from the solution of its method, or $p_{ab} = -c_{ab}$ from plates, and the
capacitance matrix of the cluster is assembled from the pairs,

$$
c_{ab} = -p_{ab} \quad (a \ne b), \qquad c_{aa} = \sum_{b \ne a} p_{ab}, \tag{108}
$$

which is symmetric with rows that sum to zero, as (71) requires. Gradient
and Hessian go through (107) and (108) unchanged: they are linear, so the
derivatives of the assembled matrix are exact. When every pair of a cluster
is dense and its mesh fits `max_elements`, the cluster is solved whole by
(70) instead. A dense window whose smallest overlap would exceed
`max_elements` is solved in windows, and its pairs are reported as such.

### Accuracy and defaults

Measured against one dense system on the same mesh, with elements of half
a gap up to 2 um and $d_i$ of 30 um, over combs of 6 fingers with overlaps
of 5 to 160 um and gaps of 2 to 8 um, single walls of 4 to 400 um across 2
to 8 um, and racks of 12 teeth 2 to 32 um wide across 1 to 4 um, as the
ratio to the dense result:

| Feature size $\lambda$ | Plates: capacitance | Plates: pull across the gap, side spring | Windows: capacitance, force |
| --- | --- | --- | --- |
| about 1 | 0.27 to 0.70 | 0.63 to 0.76 | within 1 percent |
| 2 | 0.39 to 0.72 | 0.77 to 0.89 | within 1 percent |
| 4 to 5 | 0.61 to 0.79 | 0.86 to 0.94 | within 1 percent |
| 8 to 10 | 0.75 to 0.88 | 0.93 to 0.97 | within 1.5 percent |
| 12 to 20 | 0.81 to 0.93 | 0.95 to 0.98 | within 2 percent |
| 25 to 40 | 0.88 to 0.96 | 0.98 to 0.99 | within 2 percent |
| 50 | 0.93 to 0.94 | 0.988 to 0.994 | within 2 percent |
| 80 to 100 | 0.965 to 0.975 | 0.993 to 0.997 | within 1 percent |
| 200 | 0.98 | 0.997 to 0.998 | within 1 percent |

The error of plates is the share of fringe field, and it falls with
$\lambda$ for every geometry; length alone does not order it, since a
wall hundreds of microns long across a wide trench is still mostly fringe.
At one $\lambda$ a comb is closer than a single wall, whose two open ends
carry fringe that a finger between two others does not.
Two errors of plates do not fall with $\lambda$: the force of a comb along
its fingers, which comes from the finger tips, is 0.97 to 1.11 of the field
at every overlap, and plates have no detent stiffness at all. The error of
windows grows with the number of windows and not with $\lambda$.

The defaults are $\lambda_d = 10$ and $\lambda_p = 100$. Plates take only
pairs whose walls are at least a hundred gaps long, where every geometry
measured is within 3.5 percent in capacitance and 1 percent in pull and
side spring. Below ten gaps fringe carries from an eighth to most of the
capacitance and a pair gets one dense system. Pairs between take windows,
within 2 percent on the geometries measured above, whose facing walls are
no longer than a few $d_i$; between pads several times longer than $d_i$
windows drop the far field and are 3 to 4.5 percent low, as the note
below measures. A caller trades accuracy for time with the
two cutoffs: $\lambda_p$ of 15 puts on plates combs whose fingers overlap
more than about fifteen gaps. At twenty gaps such a comb is 7 to 9 percent
low in capacitance, within 3 percent in side spring, and up to 8 percent
off in force along the fingers, found in a fraction of the time.
$\lambda_d = \lambda_p = 0$ takes every pair as plates, $\lambda_d = \lambda_p
= \infty$ solves every pair in one dense system as far as `max_elements`
allows, and $\lambda_d = 0$, $\lambda_p = \infty$ solves every pair in
windows.

!!! note "Pieces found for every cluster"
    Pieces of the plates method are found for every cluster whatever the
    methods, since they give the feature sizes. They cost a small part of
    a field solution.

!!! note "Methods fixed at rest"
    Methods and cores are fixed at rest, and windows follow travel by
    (164). Evaluated at displaced bodies, a pair keeps its method, so
    capacitance stays a smooth function of displacement, and (106) is never
    applied anew to a moved geometry.

!!! note "Reach of a window"
    A window reaches no farther than its overlap, at most $d_i$, so charge
    on parts of a conductor whose windows reach no partner is not counted.
    A cluster solved whole counts it. The two differ most between large pads
    across wide gaps: by 4.5 percent for pads 200 um long across 16 um with
    $d_i$ of 30 um, and by 3 percent with $d_i$ of 60 um.

!!! note "Element size"
    A pair solved by the field needs `element_size`; a cluster whose pairs
    are all plates is not meshed.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $\gamma(\mathbf{r})$ | gap from a point of a wall along its normal to the first conductor | `FacingWalls.gap` at the ends of each piece |
| $\Gamma_{a \to b}$ | walls of conductor $a$ that face conductor $b$ within $d_i$ | `SidewallCapacitance.walls`, `Cluster.walls` |
| $d_i$ | widest gap across which two walls count as a plate | `ElectrostaticParams.interaction_distance` |
| $m_j$, $x$ | slope of the edge faced and the argument of the piece integral | `electrostatics.facing.STEEPEST` bounds the slope; `SERIES_BELOW` the argument |
| $w_r$, $L_r$ | plate capacitance per permittivity and length of a run of pieces | internal to `electrostatics.facing` |
| $\lambda_{AB}$ | feature size of a pair of nets | `SidewallCapacitance.feature_size`, `Interaction.feature_size` |
| $\lambda_d$, $\lambda_p$ | dense and plates cutoffs | `ElectrostaticParams.dense_cutoff`, `ElectrostaticParams.plates_cutoff` |
| $m_{AB}$ | method of a pair | `SidewallCapacitance.method`, `Interaction.method` |
| $p_{ab}$ | capacitance of a pair of conductors from its method | entries of `SidewallCapacitance.capacitance` |
