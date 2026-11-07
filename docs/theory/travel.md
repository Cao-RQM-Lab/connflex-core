# Windows and rows over travel

A window chosen at rest holds the elements within its overlap, and a row the
cells that were translates of one another at rest. Bodies that travel farther
than an overlap, as the teeth of a stepper or a rack do, carry walls into and
out of both. This page gives the rule by which a window takes a new set of
elements and a row stops counting cells that travel has carried near its
ends, so that the capacitance model keeps the accuracy it had at rest over
many pitches of travel. It acts inside `CapacitanceModel.at`.

A window chosen at rest holds the elements within its overlap $\delta_k$ of
its core, (91). Bodies that travel farther than that bring walls before the
core that the window does not hold, and carry away walls it does. Teeth that
slide past each other by a few pitches, as in a stepper or a rack, would
leave each core facing open space in its window, and the pull would be tens
of percent off. Where the bodies are now, with $\mathbf{r}_i(\mathbf{q})$ the
element centers and $B_k(\mathbf{q})$ the box of the core's centers, a set
of elements $S$ still serves window $k$ while it holds every element within
$(1 - \mu)\, \delta_k$ of its core:

$$
\bigl\{ i : \mathbf{r}_i(\mathbf{q}) \in B_k(\mathbf{q}) \oplus [-(1 - \mu)\, \delta_k,\, (1 - \mu)\, \delta_k]^2 \bigr\} \subseteq S, \qquad \mu = \tfrac{1}{4}. \tag{164}
$$

A window keeps every set it has taken, the one of (91) at rest first. It
takes the set it took last while that serves, otherwise the newest that
serves, and otherwise a new one: the elements within $\delta_k$ of its core
where they are now, by (91) with the overlap found at rest, halved while the
window would hold more than `max_elements`. A window that held no pair at
rest is solved once a pair comes within $(1 - \mu)\, \delta_k$ of its core.
A set changes only after travel of $\mu\, \delta_k$ or more, and while it
holds, (92) to (95) hold as they are: $\mathbf{C}(\mathbf{q})$ is smooth and
(95) its exact derivative. Where a window changes its set,
$\mathbf{C}(\mathbf{q})$ steps by what the overlap leaves out, of the size
of the asymmetry of (94).

A row goes stale too. Its cells are elements fixed at rest, and the
periodic solution counts every element of a cell as if the row went on
without end. Travel along a row moves the elements of one body through the
cells of the other, and the runs that end the row with it. With $s_c^-$ and
$s_c^+$ where the runs of conductor $c$ of the row begin and end along it,
carried with its body to where it is now, cell $k$ is counted while every
one of its elements lies at least $d_i$ inside every run:

$$
\max_c s_c^- + d_i \;\le\; \mathbf{r}_i(\mathbf{q}) \cdot \mathbf{e} \;\le\; \min_c s_c^+ - d_i
\quad \text{for every element } i \text{ of cell } k. \tag{165}
$$

At rest every cell is counted, by the construction of (153). Cells that are
not go to windows of their own, in blocks of whole cells of at most $n_w$
elements, which take $d_i$ as their overlap and follow (164). The cell
solved for a row is the counted cell nearest its middle.

On a bar with teeth at a pitch of 10 micrometers under a fixed block with
teeth 4 micrometers from tip to tip, slid 3 to 7 pitches along with $d_i$ of
30 micrometers, windows kept from rest put the pull 17 to 41 percent off a
solution of the whole mesh, and windows and rows that follow travel 1.5 to
4.1 percent off. Windows and rows planned anew where the bar is are 1.3 to
2.3 percent off. Both drop field beyond $d_i$, which changes as the
overhang of the bar does, and with $d_i$ of 60 micrometers both stay within
1.6 percent.

!!! note "History of a window"
    Which set a window takes depends on where the bodies have been: two
    histories that reach one displacement can differ by a change of set.
    A static solution or a transient follows one path, and sees a set
    change once each time travel crosses the margin of (164).

!!! note "Cost of the check"
    Checking (164) takes one query of a tree of element centers per window.
    A surrogate's expansion, (149), checks nothing; its anchors are field
    solutions, each with the sets of its moment.

!!! note "Rows that turn"
    Rows that turn are not followed beyond (154): an arc of teeth is no row
    of translates, and windows cover it.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $\delta_k$ | overlap of window $k$, found at rest | `SidewallCapacitance.windows` counts the windows |
| $\mu$ | share of its overlap a window gives up to travel before it takes a new set | `electrostatics.windows.TRAVEL_MARGIN` |
| $d_i$ | interaction distance, the margin a counted cell keeps from the ends of a row | `ElectrostaticParams.interaction_distance` |
| $n_w$ | most elements in a window of cells a row no longer counts | `ElectrostaticParams.window_elements` |
