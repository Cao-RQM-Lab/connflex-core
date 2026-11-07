# Layout, release, and beams

The first stage turns the polygons of a layout into the parts of a mechanical
model: the layers after healing, the buried oxide that survives release, the
rigid bodies, and the beams that join them. It takes the device, handle, and
optional release polygons of a `Layout`, the layer `Stack`, and an
`ExtractionParams`, and `extract` runs it. Its results are the bodies and beams
of a `Structure`, and the candidates it rejected with the reason for each.

As a self-contained polygon-to-mechanics simulator, connflex-core automatically
extracts mechanical components given a GDS file (or a number of generic
polygons). connflex-core targets the silicon-on-insulator (SOI) process and
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
and $r_h$ is the healing tolerance. Device and handle layers are then
simplified with the Douglas–Peucker algorithm with the tolerance
`simplify_tolerance`, $r_h$ by default; a release layer is simplified without
being healed. Healing adds material and never removes it, so it leaves
protrusions and shallow notches of any width in place; simplification moves
no vertex farther than the tolerance and may remove such features up to the
tolerance deep. A drawn arc moves inward by up to the tolerance, a hundredth
of a one-micrometer gap at the default.

Release holes may be provided in a dedicated layer or within the device layer
polygon itself. In the latter case, the release holes are identified by a
single area threshold $A_\text{min}$. Holes at or below $A_\text{min}$ are
release openings: they are filled for connectivity, body assembly, and beam
extraction, enter the release estimate of (2), and are subtracted from the
device-layer mass, center of mass, and inertia of their body.

Surviving BOX regions after release, $O$, are regions protected by both device
and handle layers, subtracted by the release hole regions, then isotropically
eroded by release distance $d$:

$$
O = \bigl((D \cap H) - R\bigr) \ominus d. \tag{2}
$$

The largest handle polygon by area is selected as the anchor for the entire
structure. connflex-core iterates through all polygons in $D$ and $H$, discarding
all polygons without a connection to the anchor.

All features of width at most $w_\text{min}$ are removed through an erosion
then an expansion by $w_\text{min}/2$,

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

Beam candidates must overlap with exactly two different bodies and keep no
surviving oxide $O$ of (2) underneath in order to qualify as a real beam: a
candidate wider than $2d$ keeps oxide under its middle, is bonded to the
handle there, and is rejected with the oxide area as its reason. A discarded
candidate that overlaps exactly one body is
a stub, such as a comb finger: it is rigidly attached to that body and
contributes mass and moment of inertia, but no stiffness. Discarded
candidates that overlap no body or three or more bodies carry neither mass
nor stiffness. Once a beam is accepted, its length $L$ is calculated from the
Euclidean distance of the overlap region centroids, its width $w$ is calculated
by dividing its geometric area by its length, and its thickness $t$ is the
device layer thickness.

Only straight beams are supported. A candidate is taken as a straight strip
of uniform width between its two attachment points, and its shape is tested
against that: with less than $f_s = 0.9$ of its area inside the strip of
width $w$ along the chord, extended by $w / 2$ past either end for oblique
end faces, as a bent, forked, or stepped candidate has, such as a folded
flexure whose fold is narrower than $w_\text{min}$, it is rejected with that
reason, since taken as one straight beam along the chord it would overstate
its stiffness by a large factor. A bend is drawn as a block wider than
$w_\text{min}$, so that it is a body and each straight leg is a beam of its
own. A slanted end lies within the strip; a taper of two to one keeps 0.92
of its area inside and passes as a strip of its mean width.

!!! note "Canonical vertex order"
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

!!! note "Release holes, healing, and the oxide estimate"
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

!!! note "Width, end points, slivers, and oxide under a candidate"
    The width uses the area of the unexpanded candidate $D - B$, so the
    expansion by $\epsilon$ cannot inflate it. Beam endpoints are placed at the
    overlap centroids at the device mid-plane, $z = -t_\text{device}/2$.
    Candidate parts narrower than $r_h$ are removed by an opening before
    qualification, since the body opening in (3) does not reproduce device edges
    exactly and would otherwise leave hairlines that join separate candidates.
    Candidates whose width falls outside the configured bounds are discarded
    like any other rejected candidate. Oxide under a candidate is measured on
    the unexpanded candidate and counts when its area exceeds $\epsilon L$:
    the eroded oxide of (2) is polygonized with arcs at reflex corners, and
    where a candidate meets a body the arc's chords can cross the candidate's
    root by a sliver that is no bond. Beams carry stiffness but no mass.

!!! note "Opening polygon by polygon"
    The opening in (3) and the sliver removal on the candidates are applied
    polygon by polygon. An opening never leaves its polygon, so the parts of a
    disjoint set stay disjoint and the result is the same as opening the set as
    a whole, without the union that a single buffer over the set performs.
    Along an edge at an angle to the axes, the mitre buffers leave a vertex a
    few $10^{-15}$ m off the line wherever a feature was removed; vertices
    within $10^{-12}$ m of the line through their neighbors are dropped after
    the opening, so that an outline repeats its edges exactly where its
    features repeat and no boundary element is cut at such a vertex.

!!! note "Stubs, junctions, and mass"
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

!!! note "Moments about the center of the bounding box"
    The sums of (11b) are evaluated about the center $\mathbf{c}_0$ of the
    bounding box of the region and moved to the origin by (11c) with the shift
    reversed. About the global origin every term of (11b) is of order
    $|\mathbf{c}_0|^2 A$, whereas the central moments of (10) are of order
    $\sigma^2 A$ for a body of extent $\sigma$, so a body a few millimeters
    from the origin would lose about seven significant digits of inertia in
    the sum. Summed about $\mathbf{c}_0$, the only cancellation left is the
    shift to the center of mass in (11c).

## Symbols on this page

| Symbol | Meaning | In the code |
| --- | --- | --- |
| $D$, $H$, $R$ | device, handle, and release-opening regions | `Layout.device`, `Layout.handle`, `Layout.release`; openings in `Structure.openings` |
| $r_h$ | healing and simplification tolerance | `ExtractionParams.heal_distance`, `ExtractionParams.simplify_tolerance` |
| $A_\text{min}$ | largest area of a hole taken as a release opening | `ExtractionParams.release_hole_max_area` |
| $d$ | release undercut distance | `ExtractionParams.release_distance` |
| $O$ | surviving buried oxide | `Body.box` |
| $w_\text{min}$ | narrowest feature of a rigid body | `ExtractionParams.body_min_width` |
| $B$ | device silicon of the rigid bodies | `Body.device` |
| $P$ | beam candidates | `Structure.beams`, `Structure.rejected`, `Body.stubs` |
| $\epsilon$ | expansion of candidates to overlap their bodies | `ExtractionParams.beam_overlap_margin` |
| $L$, $w$, $t$ | length, width, and thickness of a beam | `Beam.length`, `Beam.width`, `Beam.thickness` |
| $f_s$ | share of a candidate's area in the straight strip along its chord | `geometry.beams.STRAIGHT_FILL` |
