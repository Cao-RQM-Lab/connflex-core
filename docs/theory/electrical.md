# Electrodes and links

The electrical model of a structure is a network. Its nodes are electrodes,
pieces of silicon that have one potential and one rigid motion, and its
branches are the links that conduct between them: beams, bridges, vias, and
wires. This page defines both from the polygons of the extraction and gives
the resistance of each kind of link. `extract` builds them, and they are read
from `Structure.nets`, `Structure.electrodes`, and `Structure.links`.

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
overlaps over the largest area, or with the anchor when it overlaps none.

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

!!! note "Why the node is an electrode"
    Neither a body nor a net can serve as a node. The anchor body carries every
    bonded pad and every stator, which sit on different nets, and one net runs
    through every body that beams join to its pad. An electrode is the largest
    piece of silicon that has both one potential and one rigid motion.

!!! note "A stub that bridges"
    A candidate that touches two polygons of the same body is a stub of that
    body for mass, since it moves with it, and a bridge for conduction, since
    it joins two electrodes. The silicon of a bridge or a beam belongs to no
    electrode.

!!! note "Order of the electrodes"
    Device electrodes are listed first and handle electrodes after them, each
    ordered by body and then by descending area of the body polygon.

!!! note "One resistance per via"
    A via carries one resistance rather than a value per area because the
    `vias` layer states only that device and handle are connected, not how the
    contact is made. A via that straps several device electrodes joins each of
    them to the handle through $R_v$, so two of them see $2 R_v$ between each
    other.

!!! note "No default resistivity"
    Resistivity of an SOI device layer spans orders of magnitude with doping,
    so it has no default. Without it the resistance of a beam link is left
    unknown rather than guessed, and extraction is otherwise unchanged.

!!! note "Where a wire may end"
    A wire end must lie on silicon of an electrode. An end on a beam or a
    bridge is an error, since that silicon belongs to no electrode, and so
    is an end on no conductor of its layer, which also catches silicon that
    pruning dropped. A point in a release hole of a pad counts as on the
    pad, as the pad's region is filled. A wire with both ends on one
    electrode joins nothing and is an error.

!!! note "Wires belong to the nets"
    Wires belong to the nets, so every analysis that reads nets sees bonded
    pads as one: the capacitance between nets, potentials and charges in
    static equilibrium, and the resistor network of the links.

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $B_j$, $P_k$ | device polygons of the bodies and the beam candidates | `Body.device`, `Structure.rejected` |
| $T_k$ | body polygons that candidate $k$ touches | `RejectedCandidate.bodies` |
| $E_j$ | electrode of body polygon $j$ | `Electrode.region`, `Structure.electrodes` |
| $\rho_e$ | electrical resistivity of the device layer | `Material.resistivity` |
| $R$ | resistance of a beam link | `Link.resistance` |
| $R_v$ | resistance of one via link | `ExtractionParams.via_resistance` |
