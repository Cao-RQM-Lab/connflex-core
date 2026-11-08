# Field solutions over travel

Equilibria and transients evaluate the capacitance model at many
displacements close to each other, and a field solution at each would cost as
much as the model at rest every time. This page gives two ways the field
method keeps that cost down without giving up the exactness of its
derivatives: a surrogate that reuses a kept solution within a radius around
it, and rows of identical cells that a comb repeats, solved as one cell of an
endless row. Both act inside `CapacitanceModel.at`; `surrogate_tolerance` and
`periodic_min_cells` of `ElectrostaticParams` switch them on. How windows
and rows follow bodies that travel far is on the [next page](travel.md).

## Surrogate for displaced bodies

An equilibrium or a transient evaluates the capacitance model at many
displacements close to each other, and with the field method every
evaluation solves each cluster anew at the cost of the solve at rest. Each
such solution gives more than values: with $\mathbf{m}$ the in-plane motion
of the bodies of a cluster, $(a_x, a_y, a_\psi)$ of (114) for each, it gives
the matrix $\mathbf{c}$ and its first and second derivatives in
$\mathbf{m}$ exactly, (81) and (84). With the fringe of (163) the motion
also holds the rise $u_z$ and the turns $\theta_x$, $\theta_y$ of every
body, of which the fringe and its derivatives are kept as well. Motions are
compared by how far they carry any silicon of the cluster,

$$
\rho(\mathbf{m}, \mathbf{m}') = \max_B \Bigl( \bigl\| (a_x, a_y)_B - (a_x, a_y)'_B \bigr\| + R_B\, \bigl| a_{\psi,B} - a'_{\psi,B} \bigr| + |\Delta u_{z,B}| + R_B \bigl( |\Delta \theta_{x,B}| + |\Delta \theta_{y,B}| \bigr) \Bigr), \tag{148}
$$

with $R_B$ the largest distance from the center of body $B$ to a corner of
the bounding box of any of its electrodes in the cluster, and $\Delta$ the change of
rise and turns between the two motions. A solution at
$\mathbf{m}_a$ is kept with a radius $r_a$. At a motion $\mathbf{m}$ within
the radius of a kept solution, the nearest such one supplies the cluster
by its expansion to second order, with $\boldsymbol{\delta} = \mathbf{m} -
\mathbf{m}_a$:

$$
\mathbf{c}(\mathbf{m}) \approx \mathbf{c}_a + \frac{\partial \mathbf{c}_a}{\partial \mathbf{m}}\, \boldsymbol{\delta} + \tfrac{1}{2}\, \boldsymbol{\delta}^\top \frac{\partial^2 \mathbf{c}_a}{\partial \mathbf{m}^2}\, \boldsymbol{\delta}, \qquad
\frac{\partial \mathbf{c}}{\partial \mathbf{m}}(\mathbf{m}) \approx \frac{\partial \mathbf{c}_a}{\partial \mathbf{m}} + \frac{\partial^2 \mathbf{c}_a}{\partial \mathbf{m}^2}\, \boldsymbol{\delta}, \qquad
\frac{\partial^2 \mathbf{c}}{\partial \mathbf{m}^2}(\mathbf{m}) \approx \frac{\partial^2 \mathbf{c}_a}{\partial \mathbf{m}^2}. \tag{149}
$$

Elsewhere the cluster is solved anew and the solution kept. The gradient is
the exact derivative of the expanded matrix, so forces within one radius
derive from one energy. The error of the gradient grows as $\rho^2$, which
sets the radii: when a new solution at $\mathbf{m}_b$ lies within four
radii of the nearest kept one, the error $e$ of that one's expansion there,
the largest of the relative errors of $\mathbf{c}$ and of its gradient in
the Frobenius norm, and of the fringe matrix and its gradient when the
fringe is kept, resizes it, and the new solution starts with the
result:

$$
r_a \leftarrow \min\Bigl( 2\, r_a,\; 0.5\, \rho(\mathbf{m}_b, \mathbf{m}_a) \sqrt{\tau_c / e} \Bigr), \qquad
r_b = r_a, \qquad
r_\text{first} = \sqrt{\tau_c}\; g_\text{min}, \tag{150}
$$

with $\tau_c$ the tolerance and $g_\text{min}$ the narrowest gap between
electrodes of different bodies in the cluster at rest. A new solution
farther than four radii from the nearest kept one inherits that one's
radius unchanged. A radius at most
doubles per check, so it grows only as fast as evidence supports it, and
a check along a step without motion in the plane can shrink a radius but
never grows it: rise and tilt leave the sidewalls unchanged and say nothing
about their expansion in the plane. A solution at a motion where one is
already kept, as an exact evaluation or an equilibrium solved anew makes,
replaces the kept one and keeps its radius.

!!! note "The tolerance is a target"
    The tolerance is a target, not a bound. The radius rule assumes the
    error of an expansion grows as $\rho^2$ in every direction from the
    one it measured, and the factor of one half leaves room for directions
    in which it grows faster; a radius is never tested against an exact
    solution before it is used. On a comb that closes its side gaps by a
    fifth the gradient stays within 0.6 of the tolerance at $\tau_c =
    10^{-2}$ and within a quarter of it at $10^{-3}$; with the factor at
    0.8 it reached 1.4 times the tolerance. A caller who needs a bound
    evaluates with `exact=True` where it matters.

!!! note "Second derivative held over a radius"
    The second derivative is held constant over a radius, so its error is
    of the order of $\sqrt{\tau_c}$, not $\tau_c$. It enters only Newton and
    integrator Jacobians, which change how fast a solve converges and not
    where it ends. An equilibrium therefore solves its result anew, so its
    stiffness and model are exact at its displacement, which converged on
    the expanded forces and carries an error of order $\tau_c$; a transient
    keeps the expanded values. Because the second derivative is held, its
    map onto the degrees of freedom of the bodies is kept with the solution
    and shared by every evaluation that expands it.

!!! note "Steps between kept solutions"
    Moving from one kept solution to the next changes the capacitance by up
    to about $\tau_c$ times itself. A transient's energy ledger counts that
    change, so its imbalance reflects the tolerance of the surrogate as well
    as that of the integrator.

!!! note "Kept solutions belong to the model"
    Kept solutions belong to the capacitance model at rest and every model
    evaluated from it, so a sequence of equilibria and transients built on
    one model reuses them all. Which solution supplies a motion depends on
    which were kept before, within the tolerance.

!!! note "Plates evaluated exactly"
    Clusters whose pairs are all taken as plates are cheap and evaluated
    exactly; so are plates between device and handle.

## Rows of identical cells

Combs repeat one motif many times, and a field problem solved in windows
then solves nearly the same window over and over. A row repeated without
end along a period $p\, \mathbf{e}$, with $\mathbf{e}$ a unit vector, needs
one cell: summed over all copies, the kernel of (61) becomes

$$
G_p(w) = -\frac{1}{2\pi} \ln \Bigl| 2 \sin \frac{\pi w}{p} \Bigr|, \qquad
w = \bar{e}\, (z - z'), \tag{151}
$$

with $z$ and $z'$ points as complex numbers and $\bar{e}$ the conjugate of
$\mathbf{e}$ as one, which turns the period onto the real axis. It differs
from the sum of (61) over the copies by a constant, which a row with no net
charge per cell does not feel. The neutrality of (63) holds per cell. For
the element integrals, $G_p$ splits into the nearest copies and a smooth
remainder,

$$
\ln \Bigl| 2 \sin \frac{\pi w}{p} \Bigr| = \text{const} + \sum_{|k| \le K} \ln |w - k p| + R_K(w), \qquad
R_K(w) = \ln \Bigl| \frac{\sin(\pi w / p)}{\pi w / p} \Bigr| - \sum_{k=1}^{K} \ln \Bigl| 1 - \frac{w^2}{k^2 p^2} \Bigr|, \tag{152}
$$

with $K = 2$. The second element of each pair is first taken to its copy
nearest the first, which (151) allows, so that $|\operatorname{Re} w|$
stays within a period or so wherever the bodies have carried the elements of
a cell. The copies enter through the closed form (67) and its derivatives,
and $R_K$, analytic for $|\operatorname{Re} w| < (K+1)\,p$, by four-point
Gauss-Legendre quadrature on each element of a pair, with its
first and second derivatives from those of $\ln \sin$. Every copy of an
element moves with it and turns about the copy of its pivot, so the pairs
of elements of one body enter its derivatives in turn, which (84) leaves
out: the copies of a body turn about different points and move against
each other, while under travel they move together and those terms cancel.
Where a wall crosses the cell boundary, the elements on either side of it
are copies that slide against each other under such a turn, so the second
derivative of one cell in turn is logarithmically singular there; the
closed form gives a finite value, which enters only the second-order term
of (154), since the Hessian in turn comes from the window.
The matrix of (64) for one cell is symmetric but not always definite, and
is factored with pivoting.

A row is found in the outlines. Where an outline repeats its edges,
translated by a vector $\mathbf{t}$, over at least $n_p$ motifs, the corners
of that run belong to a row of period $\mathbf{t}$, with the shortest motif
kept where runs of several motifs cover the same edges. Runs of two or more
conductors with one period whose extents along it overlap, and lie within
$d_i$ of each other across it, form one row: along $\mathbf{t}$ it spans
where every one of those conductors repeats, less $d_i$ at either end, and
across it the band of their runs. Cells are slabs of width $|\mathbf{t}|$
across that band, whose boundaries pass through the widest gap between
element midpoints modulo $|\mathbf{t}|$, and an element belongs to the cell
that holds its midpoint. The interior of the row is the longest run of
cells, at least $n_p$, whose elements are translates of those of the first,

$$
\mathbf{r}_{e + m} = \mathbf{r}_e + k\, \mathbf{t} \quad \text{for every element } e \text{ of cell } 0 \text{ and its counterpart } e + m \text{ in cell } k, \tag{153}
$$

to a thousandth of the shortest element. Its elements leave the cores of
the windows but stay in the mesh, so windows near the row see its field.

The middle cell $m$ of the row is solved with (151) at the bodies where they
are. It gives the capacitance $\mathbf{c}$ that one cell holds and its
derivatives $\mathbf{g}$ and $\mathbf{H}$ in the travel $\boldsymbol{\tau}_B$
and turn $\psi_B$ of each body about the cell center. With a body turned by
$\psi_B$, cell $k$ sees that body shifted from where the middle cell sees it
by $\boldsymbol{\delta}_{kB} = \bigl( \mathbf{R}(\psi_B) - \mathbf{I} \bigr)
(\mathbf{c}_k - \mathbf{c}_m)$, with $\mathbf{c}_k$ the cell centers at
rest, and holds

$$
\mathbf{c}_k \approx \mathbf{c} + \mathbf{g}\, \boldsymbol{\delta}_k + \tfrac{1}{2}\, \boldsymbol{\delta}_k^\top \mathbf{H}\, \boldsymbol{\delta}_k + \boldsymbol{\psi}^\top \bigl( \Delta\mathbf{g} + \Delta\mathbf{H}\, \boldsymbol{\delta}_k \bigr) - \tfrac{1}{2}\, \boldsymbol{\psi}^\top \Delta\mathbf{H}\, \boldsymbol{\psi}, \qquad
\frac{\partial \mathbf{c}_k}{\partial \boldsymbol{\ell}} \approx \mathbf{g} + \mathbf{H}\, \boldsymbol{\delta}_k + \mathbf{P}_\tau\, \Delta\mathbf{H}\, \boldsymbol{\psi}, \tag{154}
$$

with $\boldsymbol{\psi}$ the turns of the bodies from rest on the entries
of $\boldsymbol{\ell}$ in turn and zero on those in travel,
$\Delta\mathbf{g}$ and $\Delta\mathbf{H}$ the excess of the window over
the periodic solution described below, and $\mathbf{P}_\tau$ the
projection on the entries in travel.

The row holds the sum over its cells. A rigid motion of body $B$ about its
pivot $\mathbf{P}_B$, travel $\mathbf{u}_B$ then turn $\psi_B$, moves the
center of cell $k$ by $\mathbf{u}_B + (\mathbf{R}(\psi_B) - \mathbf{I})
(\mathbf{c}'_k - \mathbf{P}_B)$ and turns it by $\psi_B$, with
$\mathbf{c}'_k$ where that center is now. So with $\mathbf{a}_{kB} =
\mathbf{c}'_k - \mathbf{P}_B$ and the local parameters
$\boldsymbol{\ell} = (\boldsymbol{\tau}_B, \psi_B)$ of every body,

$$
\frac{\partial \boldsymbol{\tau}_{kB}}{\partial \mathbf{u}_B} = \mathbf{I}, \quad
\frac{\partial \boldsymbol{\tau}_{kB}}{\partial \psi_B} = \hat{\mathbf{z}} \times \mathbf{a}_{kB}, \quad
\frac{\partial^2 \boldsymbol{\tau}_{kB}}{\partial \psi_B^2} = -\mathbf{a}_{kB}, \qquad
\nabla^2 \mathbf{c}_k = \mathbf{J}_k^\top \mathbf{H}\, \mathbf{J}_k - \sum_B \frac{\partial \mathbf{c}_k}{\partial \boldsymbol{\tau}_B} \cdot \mathbf{a}_{kB}\; \mathbf{e}_{\psi_B} \mathbf{e}_{\psi_B}^\top, \tag{155}
$$

with $\mathbf{J}_k$ the Jacobian of the local parameters in the body motions.

One cell cannot carry the derivatives in turn. A body that turns rigidly
about a cell center also shifts the copies of that cell across the row, each
by its distance times the turn; no motion of a single periodic cell does
that. For the derivatives in $\psi_B$, $\mathbf{g}$ and $\mathbf{H}$
therefore take their entries from a window around the middle cell, of
$d_i$ or two periods, whichever is larger, which counts charge on that
cell and turns its bodies rigidly about the cell center. The periodic
solution keeps the capacitance and the derivatives in travel. Its
capacitance misses those same shifts: solved at the bodies where they are,
it has every copy turned about its own center, so with the bodies turned
by $\boldsymbol{\psi}$ from rest it is short of the capacitance whose
derivative the window gives by the integral over the turn of the excess
$\Delta\mathbf{g}$ of the window's $\mathbf{g}$ over the periodic one on
the entries in turn. The last two terms of (154) restore that integral to
second order, with $\Delta\mathbf{H}$ the same excess in $\mathbf{H}$ on
the entries with an index in turn, so that the derivative of $\mathbf{c}_k$
in $\boldsymbol{\psi}$ is $\mathbf{g} + \mathbf{H}\, \boldsymbol{\delta}_k$
to second order in $\boldsymbol{\psi}$ and the torque is the derivative of
the energy; the term $\mathbf{P}_\tau \Delta\mathbf{H}\, \boldsymbol{\psi}$
is the derivative in travel of $\boldsymbol{\psi}^\top \Delta\mathbf{g}$.
On a comb of 40 cells turned by $10^{-4}$ the torque was 7.7 percent off
the derivative of the capacitance without them, and the capacitance
changed with the turn by as much against one dense solution.

What a row leaves to the windows repeats with it: congruent rows have
congruent ends, and the walls beside their cells repeat with the cells.
Cut by medians over all of a mesh, these elements fall into cores that
differ from row to row, so they are cut in the frame of their row instead.
With $s = \mathbf{r} \cdot \mathbf{e}$ along a row, $s_-$ where its first
cell begins and $s_+$ where its last ends, and $\ell = \lceil d_i / |\mathbf{t}|
\rceil\, |\mathbf{t}|$ the fewest whole cells that span $d_i$, an element
within $2 d_i + |\mathbf{t}|$ of the cells of a row joins the nearest such
row, in group

$$
z = \begin{cases}
\text{before} & s < s_-, \\
\text{beyond} & s > s_+, \\
\bigl\lfloor (s - s_-) / \ell \bigr\rfloor & \text{otherwise},
\end{cases} \tag{156}
$$

and each group is cut into cores by medians, as in (91). Elements near no
row form one more group. Rows that are translates of one another, with
surroundings that are too, then give cores and windows that are translates,
which are solved once.

!!! note "Ends of a real row"
    A real row ends, and its cells near either end differ from those of an
    endless row. The interior stops $d_i$ short of where the row ends, and
    the windows of (156) take what lies beyond. Against one dense solution
    of a comb of 40 cells, at rest and with the shuttle moved and turned,
    rows with $d_i$ of 10 micrometers put capacitance within $6 \times
    10^{-5}$, the pull along the fingers within $1 \times 10^{-3}$, the
    spring in turn within $8 \times 10^{-4}$, and the side spring within
    $2 \times 10^{-5}$. Windows with an overlap of 30 micrometers are off by
    $3 \times 10^{-4}$, $4 \times 10^{-3}$, $3 \times 10^{-3}$, and $3
    \times 10^{-5}$ and take twenty times as long; with 10 micrometers, the
    pull is off by $1.5 \times 10^{-2}$. On a comb of 80 cells rows put
    capacitance $5 \times 10^{-4}$ high and windows of 30 micrometers $1.5
    \times 10^{-3}$.

!!! note "Period of a row"
    The period of a row is the translation measured on the runs of corners
    that repeat, as drawn, taken to the grid of $10^{-10}$ m on which runs
    are grouped only when it lies within rounding of a grid point. A period
    rounded to that grid would miss the cells of a row at an angle to the
    axes, or at a pitch off the grid, by a growing amount, since the translate
    test below holds the elements of every cell to a thousandth of the
    shortest element.

!!! note "Matching the elements of a cell"
    The elements of a cell are matched to those of the first cell in one
    order: by conductor, then by position along and across the row, each
    rounded to the tolerance of the translate test, a thousandth of the
    shortest element. Unrounded positions would let two elements that lie
    at one position along the row up to rounding, such as the tips of
    fingers on either side of a spine, change places from one cell to the
    next and fail the test.

!!! note "Remainder far across the row"
    Across the row, $R_K$ and its derivatives are evaluated from $\ln \sin$
    directly while $|\operatorname{Im}(\pi w / p)| \le 20$, and from the
    asymptotic forms $\ln |\sin x| = |\operatorname{Im} x| - \ln 2$,
    $\cot x = -i \operatorname{sign}(\operatorname{Im} x)$, and
    $1 / \sin^2 x = 0$ beyond, which hold there to rounding; the sine itself
    would overflow a few hundred periods across, where a body carried far
    out of its row would otherwise return no finite capacitance.

!!! note "Where rows save time"
    Rows cost least where windows would solve many copies of one window:
    a row is one cell, and translated windows are one solution. With
    bodies turned, windows that were translates at rest no longer are and
    are each solved, and rows keep the gain of their cells.

!!! note "Cells that are not translates"
    The cell mesh is the mesh of the outlines itself: an outline meshes
    each edge from its length and the angles at its ends, so congruent
    edges carry congruent elements. A row whose cells are not translates of
    each other, for instance one along a long straight edge that its
    elements do not divide evenly, is left to the windows. So is an outline
    whose repeating unit has more than 16 edges, as fingers with rounded
    tips drawn as many segments: runs are sought over motifs of up to 16
    edges, compared on a grid of $10^{-10}$ m.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $\mathbf{m}$, $\rho(\mathbf{m}, \mathbf{m}')$ | motion of the bodies of a cluster and the distance between two motions | internal to `electrostatics.capacitance` |
| $r_a$, $\tau_c$ | radius of a kept solution and the tolerance it aims at | `Cluster.anchors` counts kept solutions; `ElectrostaticParams.surrogate_tolerance` |
| $0.5$ in (150) | factor on the radius | `electrostatics.capacitance.SAFETY` |
| $G_p$, $R_K$, $K$ | periodic kernel, its remainder, and the copies taken exactly | `electrostatics.periodic`, `electrostatics.periodic.IMAGES` |
| $\mathbf{t}$, $n_p$ | period of a row and the fewest cells that make one | `ElectrostaticParams.periodic_min_cells`; rows counted in `Cluster.windows` |
| $\mathbf{g}$, $\mathbf{H}$, $\Delta\mathbf{g}$, $\Delta\mathbf{H}$ | derivatives of one cell and the excess of the window over the periodic cell | internal to `electrostatics.periodic.solve_rows` |
