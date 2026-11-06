# Usage

All quantities are SI. GDS, KLayout, and gdsfactory readers apply the unit
stored in the file or database; the shapely reader takes a `unit` scale
factor, for example `unit=1e-6` for polygons drawn in micrometers.

## Load a layout

A layout has four polygon roles: `device`, `handle`, and optionally `release`
and `vias`. Every polygon marks where the named feature is present; see the
[input conventions](conventions.md).

=== "GDS file"

    ```python
    import connflex_core as cf

    layout = cf.load_gds(
        "resonator.gds",
        layers=cf.LayerMap(device=(1, 0), handle=(2, 0), release=(3, 0), vias=(4, 0)),
    )
    ```

    File units are applied. Pass `cell="NAME"` when a file has more than one
    top-level cell. Holes stored as cut lines are repaired into polygons with
    holes, and each repaired polygon must keep its area as recorded in the
    file; a mismatch is an error naming the polygon.

=== "KLayout"

    ```python
    import klayout.db as kdb

    import connflex_core as cf

    source = kdb.Layout()
    source.read("resonator.gds")
    layout = cf.from_klayout(source, layers=cf.LayerMap(device=(1, 0), handle=(2, 0)))
    ```

    Takes a `Layout` or a `Cell`. Database units are applied and instances are
    flattened. Pass `cell="NAME"` when a layout has more than one top cell.

=== "gdsfactory"

    ```python
    import gdsfactory as gf

    import connflex_core as cf

    component = gf.import_gds("resonator.gds")
    layout = cf.from_gdsfactory(component, layers=cf.LayerMap(device=(1, 0), handle=(2, 0)))
    ```

    Database units are applied and instances are flattened.

=== "shapely"

    ```python
    import connflex_core as cf
    from shapely import box

    layout = cf.Layout.from_shapely(
        device=[box(0, 0, 200, 50), box(90, 50, 110, 150), box(0, 150, 200, 200)],
        handle=[box(-50, -50, 250, 250)],
        unit=1e-6,
    )
    ```

Polygons within one role are merged where they overlap or touch, invalid
polygons are rejected with role and index named, and `layout.device` and friends are shapely
`MultiPolygon` objects in meters. Each polygon is stored in canonical vertex order, with
the exterior ring clockwise from its smallest vertex and holes counterclockwise, so a
layout drawn with different ring start points, orientations, or hole order extracts to
bit-identical results.

Readers are looked up by name through `cf.adapters.get(name)`; `cf.adapters.available()`
lists them. A package adds a reader with `cf.adapters.register(name, function)` or by
publishing an entry point in the `connflex_core.adapters` group.

## Describe the process

A stack lists the three SOI layers from top down. Thickness is required for
each; material presets can be overridden field by field. Beam stiffness uses
the device-layer material; other layers contribute mass only.

```python
stack = cf.Stack(
    device=cf.Layer(thickness=25e-6, material=cf.Material.silicon(youngs_modulus=130e9)),
    box=cf.Layer(thickness=2e-6, material=cf.Material.oxide()),
    handle=cf.Layer(thickness=400e-6, material=cf.Material.silicon()),
)
```

## Choose parameters

Three extraction values have no default because they are process and design
specific. Everything else is listed on the [parameters page](parameters.md).
Invalid values raise `ValueError` naming the field.

```python
params = cf.ExtractionParams(
    release_hole_max_area=100e-12,
    release_distance=5e-6,
    body_min_width=10e-6,
)
options = cf.SolveOptions(n_modes=10)
```

Every parameter object is a plain frozen dataclass built from keyword
arguments, and each function takes the objects it needs as arguments. The
library reads no settings files: an application that keeps its settings in a
file parses that file itself and passes each table on, for example
`cf.ExtractionParams(**table)`.

## Extract the structure

```python
structure = cf.extract(layout, stack, params)
structure = cf.extract(layout, stack, params, wires=bonds)  # with connections off the die

structure.bodies  # rigid bodies, anchor first; polygons per layer in meters
structure.bodies[1].inertial  # mass, center of mass, inertia tensor, 6 by 6 mass matrix
structure.bodies[1].stubs  # device silicon attached to this body alone; mass, no stiffness
structure.beams  # beams joining exactly two bodies: length, width, thickness, end points
structure.rejected  # beam candidates that did not qualify, each with a reason
structure.graph  # networkx MultiGraph, bodies as nodes and beams as edges
structure.nets  # electrical nets of remaining device and handle polygons, wires included
structure.electrodes  # silicon with one potential and one rigid motion: nodes of the electrical model
structure.links  # beams, bridges, vias, and wires that conduct between electrodes
structure.openings  # release openings: filled release holes and release-layer polygons
structure.stack, structure.params  # layer stack and parameters of the extraction
structure.to_dict()  # JSON-serializable tree with GeoJSON-style polygons
```

Stages follow the [theory of operation](theory.md): release holes are filled,
device and handle layers are healed and simplified, cavities that healing
closes are filled as well, the holes and the `release` layer give the
surviving buried oxide, an anchor handle polygon is selected, polygons with no connection to it through the stack are dropped,
bodies are assembled by overlap through the opened device layer, oxide, and
handle, and device regions removed by the opening become beams when they
touch exactly two bodies and their width lies within `beam_width_min` and
`beam_width_max`. Each body carries its mass, center of mass, and inertia
tensor, computed layer by layer from the stack with exact closed-form edge
integrals and with release openings removed from the device layer;
`body.inertial.matrix` is the 6 by 6 mass matrix in body degrees of freedom.

A candidate that touches one body is a stub, one that touches three or more
is a junction; both appear in `structure.rejected` with that reason and the
bodies touched, which is the first place to look when an expected beam is
missing. A stub, such as a comb finger, is also attached to its body: it moves
with it and adds mass and inertia but no stiffness, and it is listed in
`body.stubs` rather than in `body.device`. A junction adds nothing.

Polygons in `body.device`, `body.stubs`, and `electrode.region` are filled:
release openings are modeled as perforations and are kept apart in
`structure.openings`, one disjoint polygon per opening. Mass properties
already exclude them. The JSON tree does not list them.

Mechanical analysis ([solve](#solve)) and electrostatic analysis
([capacitance model](#capacitance-model)) both work from this structure, and
each runs without the other: neither needs the other's parameters or
results. [Sidewall capacitance](#sidewall-capacitance) takes polygons
directly and needs no structure. [Static equilibrium](#static-equilibrium)
under voltage is electromechanical and takes both.

## Electrical nets

```python
nets = structure.nets
nets.count  # number of nets
nets.conductors  # connected device regions, then handle regions, each with its net
nets.vias  # each via with the conductors it joins on both sides
nets.wires  # each wire with the conductors under its two ends
nets.graph  # networkx MultiGraph, conductors as nodes, vias and wires as edges
nets.conductor_at((x, y), "device")  # conductor covering a point

structure.beams[0].net  # net of the beam's device conductor
structure.bodies[0].nets  # net of each device polygon of the body
```

Each connected device region and each connected handle region is one
conductor. A via joins every conductor it overlaps on both sides, and nets
are the connected components. A via that overlaps no conductor on one side
is ignored with a warning. Polygons dropped by anchor pruning are not
conductors. `cf.find_nets` computes the same from any lists of polygons.

### Connections off the die

```python
wires = [
    cf.Wire((1200e-6, 400e-6), (3100e-6, 400e-6)),  # bond between two pads
    cf.Wire((1200e-6, 400e-6), (0.0, 0.0, "handle"), resistance=0.5),  # pad to handle
]
structure = cf.extract(layout, stack, params, wires=wires)
```

Pads that are apart in silicon are often joined outside the device: by wire
bonds between pads or to a package lead, by traces on a board, or by the
package tying the handle to a pad. `wires` states those connections. Each
end is a point in meters on the pad, on the device layer unless a layer is
given as a third entry; any point on the pad's silicon will do, also one in
a release hole. `resistance` is that of the bond or trace and its contacts,
zero unless given.

Nets include the wires, so bonded pads are one net everywhere: in
`structure.nets`, in `electrode.net`, in `model.net_matrix`, and in the
potentials and charges of `structure.equilibrium`, where one potential on
the net drives every pad on it. Each wire is also a link of kind `"wire"`
with its resistance. An end that lies on no silicon of its layer, on a beam
or a bridge, or a wire with both ends on one electrode raises an error that
names the wire. Changing the wiring means extracting again. Many wires are
built in a loop from whatever list an application keeps:
`[cf.Wire(a, b) for a, b in bonds]`.

## Electrodes and links

```python
structure.electrodes  # device electrodes, then handle electrodes
structure.links  # beams, bridges, vias, and wires that conduct between electrodes
structure.electrode_at((x, y), "device")  # electrode covering a point

electrode = structure.electrodes[0]
electrode.body, electrode.conductor, electrode.net
electrode.region  # body polygon followed by its stubs, in meters

link = structure.links[0]
link.electrodes  # indices of the two electrodes joined
link.kind  # "beam", "bridge", "via", or "wire"
link.source  # index in structure.beams, structure.rejected, structure.nets.vias, or .wires
link.resistance  # ohms; None for a beam without resistivity
```

An electrode is silicon of one conductor that moves with one body: a device
polygon of a body with the stubs that touch it alone, or a handle polygon. It
is the node of the electrical model, since the anchor body spans several nets
and one net spans several bodies. A point on a comb finger finds the electrode
of the pad or shuttle that finger is attached to.

Links follow from the structure. A beam links the electrodes it joins. A
rejected candidate that touches several body polygons, such as a narrow strip
between two bonded pads or a junction, is a bridge between their electrodes.
A via links the device and handle electrodes it overlaps, and a wire the
electrodes under its two ends. Electrodes joined through links always share
a net.

A beam link has resistance $\rho_e L / (w t)$ from the beam's own length,
width, and thickness once the device material carries a resistivity, for
example `cf.Material.silicon(resistivity=1e-4)`. Resistivity has no default
because it depends on doping. A via link has `via_resistance` from the
extraction parameters: one value for every via, whatever its area, and
0 $\Omega$ unless set, for example `cf.ExtractionParams(..., via_resistance=2.0)`.
Bridges are ideal, and resistance inside an electrode is neglected. Exact
definitions are in
[electrodes and links](theory.md#electrodes-and-links).

## Capacitance model

```python
params = cf.ElectrostaticParams(element_size=2e-6, interaction_distance=30e-6)
model = structure.capacitance(params)

fast = cf.ElectrostaticParams(element_size=2e-6, interaction_distance=30e-6, plates_cutoff=15.0)

model.matrix  # n by n over structure.electrodes, farads
model.oxide, model.released, model.sidewall  # the three parts it sums
-model.matrix[a, b]  # capacitance between electrodes a and b
model.net_matrix  # between nets, electrodes of each net tied together
model.clusters  # groups of device electrodes whose sidewalls were taken together
model.clusters[0].interactions  # each pair of nets: feature size and method that found it
model.plates  # plates behind the oxide and released parts

model.dofs  # (body, "ux" | "uy" | "uz" | "tx" | "ty" | "tz") per degree of freedom
model.gradient  # n by n by len(dofs): change of the matrix with each
force = model.force(potentials)  # one potential per electrode
softening = model.stiffness(potentials)  # adds to the mechanical stiffness

moved = model.at(displacement)  # the same model with bodies displaced from rest
moved.matrix, moved.gradient, moved.force(potentials)  # all at that displacement
moved.clearance, moved.closest  # narrowest gap against its width at rest, and where

model.to_dict()  # JSON-serializable report: pairs, nets, clusters, what was left out
```

`structure.capacitance` collects every capacitance between electrodes at
rest into one matrix: plates between device and handle electrodes, through
oxide and across released gaps, and sidewalls between device electrodes. It
needs no mechanical solution. Rows sum to zero, so only differences of
potential matter and nothing outside the structure holds charge.

Sidewalls are found in clusters. Device electrodes whose silicon lies
within `interaction_distance` of each other, directly or through a chain of
neighbors, form a cluster; electrodes in different clusters share no
sidewall capacitance.

Within a cluster each pair of nets whose electrodes lie within
`interaction_distance` takes its own method, chosen by the feature size of
their interaction: the length of the walls that face each other in units of
the gap between them, weighted by where the capacitance is. A comb has the
finger overlap over the finger gap, a plate actuator its length over its
gap, a toothed motor its tooth width over the gap at the tips.

| Feature size | Method | Accuracy against one dense system |
| --- | --- | --- |
| below `dense_cutoff`, 10 | field in one dense system | reference |
| up to `plates_cutoff`, 100 | field in windows | within about 2 percent |
| from `plates_cutoff` | facing walls as parallel plates | within about 3.5 percent in capacitance and 1 percent in pull and side spring at 100 gaps, better beyond |

Plates cost a small part of a field solution and leave out all fringe
field, which is why they take only long walls by default. Lower
`plates_cutoff` to put more pairs on plates when time matters more than a
fringe of 5 to 15 percent; 0 for both cutoffs puts every pair on plates,
infinity for both solves every pair in one dense system, and
`dense_cutoff=0` with `plates_cutoff=float("inf")` every pair in windows.
Pairs asked to be dense whose problem needs more than `max_elements`
boundary elements are solved in windows instead: cores of at most
`window_elements` elements, each solved together with the elements around
it, out to an overlap that doubles until the charge on the core settles
within `window_tolerance` or the overlap reaches `interaction_distance`.
Measured accuracy is in [accuracy and defaults](theory.md#accuracy-and-defaults)
and what plates miss in [what plates leave out](theory.md#what-plates-leave-out).

Each cluster reports `interactions`, the feature size and method of every
pair of its nets, `windows` and `asymmetry`, which is rounding for one
dense system and a measure of the error for windows, and `walls`, the
pieces of wall that face another electrode, for a check against the
layout. `element_size` is needed as soon as one pair is solved by the
field; a structure whose pairs are all plates needs none.

The model also holds how the matrix changes when bodies move, in the six
degrees of freedom of every movable body and in the order of `result.dofs`
from `modal`. `force` returns the generalized forces for given electrode
potentials and `stiffness` the electrostatic stiffness, whose negative
diagonal entries soften the suspension. Potentials are per electrode;
`np.asarray(per_net)[model.net]` spreads one value per net. Sidewall forces
act at mid-height of the device layer, so they also tilt a body whose center
of mass lies lower.

`model.at(displacement)` evaluates the model with the bodies displaced: one
value per degree of freedom in `model.dofs`, in meters and radians from
rest. It returns a model of the same kind whose matrices, gradient, force,
and stiffness are those at the displaced bodies, so forces away from rest
are exact derivatives of the capacitance away from rest. Plates are
re-evaluated in closed form and sidewalls are found again on the moved
outlines, with the clusters, methods, mesh, and windows of the model at
rest; one call costs about as much as the sidewalls of the model at rest. Bodies that
touch after the motion raise an error.

`to_dict` gives a tree for review. It lists every pair of electrodes with
capacitance between them, split into oxide, released, and sidewall parts,
with their nets and the change of the total with each body that moves it;
the matrix between nets; the clusters with their size, windows,
asymmetry, and the feature size and method of each pair of nets; and
`unassigned_area`, the device silicon of beams and bridges,
which belongs to no electrode and is not in the model.

Either part can be left out, and then needs none of its parameters:
`cf.ElectrostaticParams(include_sidewalls=False)` gives plates to the handle
alone with no element size or interaction distance, and
`include_plates=False` gives sidewalls alone. Definitions
are in [capacitance of a structure](theory.md#capacitance-of-a-structure).

## Capacitance to the handle

```python
plates = structure.plates()  # one per device electrode over a handle electrode
plates = structure.plates(cf.ElectrostaticParams(ambient_permittivity=1.0006))

plate = plates[0]
plate.electrodes  # (device electrode index, handle electrode index)
plate.bodies  # bodies those electrodes move with
plate.oxide  # farads through surviving buried oxide; fixed
plate.released  # farads across the released gap at rest

value, gradient, hessian = plate.evaluate()  # at rest
value, gradient, hessian = plate.evaluate(device=[0, 0, -0.5e-6, 0, 0, 0])
```

A plate is the overlap of one device electrode with one handle electrode. The
part over surviving oxide is a fixed capacitor with the dielectric constant
of the box material, `cf.Material.oxide()` carrying 3.9. The rest is a gap as
high as the box layer is thick, filled with a medium of
`ambient_permittivity`, and release openings are left out of its area. A
bonded pad has both parts: oxide under its middle and the undercut rim around
it.

`evaluate` takes the six degrees of freedom of the body of each electrode,
ordered ux, uy, uz, tx, ty, tz, and returns capacitance with its gradient and
Hessian over those twelve, device body first. Vertical travel is exact,
tilt is carried to second order, and in-plane travel changes the overlap to
first order wherever the outline of the handle passes under the plate. With a
potential difference `V` across the plate, `0.5 * V**2 * gradient` is the
generalized force on the two bodies and `-0.5 * V**2 * hessian` the stiffness
it adds. A plate whose electrodes move with one body cannot change:
evaluated with the same motion for both, it returns its rest value. Beams and
bridges belong to no electrode and are not counted. Definitions and limits
are in [plates between device and handle](theory.md#plates-between-device-and-handle).

## Sidewall capacitance

```python
params = cf.ElectrostaticParams(element_size=1e-6, interaction_distance=30e-6)
result = cf.sidewall_capacitance([finger, slot], params)

result.capacitance  # n by n, farads per meter of thickness
-result.capacitance[0, 1] * 25e-6  # farads between the two across a 25 um layer
result.gradient  # n by n by bodies by 3: change with x, y, and rotation of each body
result.hessian  # n by n by bodies by 3 by bodies by 3: second derivatives
result.centers  # point each body rotates about
result.feature_size  # n by n: facing wall length over gap for each interacting pair
result.method  # n by n: "dense", "windowed", or "plates" for each interacting pair
result.charge  # charge on each element per volt on each conductor, where counted
result.mesh.start, result.mesh.end, result.mesh.conductor  # boundary elements
result.windows, result.asymmetry  # field solutions it took, and their mismatch
result.walls.start, result.walls.end  # pieces of wall facing another conductor
result.walls.conductor, result.walls.facing, result.walls.gap

result = cf.sidewall_capacitance(conductors, params, bodies=[0, 0, 1], centers=centers)
```

`sidewall_capacitance` finds capacitance between conductors drawn in the
layout plane, each pair within `interaction_distance` by the method its
feature size selects, as for the capacitance model with every conductor
its own net. Each conductor is a shapely polygon or multipolygon in meters and is one equipotential;
conductors must not overlap or touch. The problem is two-dimensional, which
holds where the layer is thick against its trenches, so the result is per
meter of thickness. Only differences of potential matter: every row of
`capacitance` sums to zero, and `-capacitance[a, b]` is the capacitance
between conductors `a` and `b`.

A field solution carries no net charge and lets the potential far away
float. As plates a wall counts toward two conductors where its normal first
meets the other within `interaction_distance`. `walls` lists the pieces of
wall that face another conductor, each with the conductor it faces and the
gap at its two ends; they give the feature sizes and the plates. Going from
`start` to `end` of a piece, the wall it faces lies to the right, so the
four points `start`, `end`, and each moved to the right by its gap outline
the plate. `charge` and `mesh` are empty when every pair is plates.

`gradient[a, b, k]` holds the derivatives of `capacitance[a, b]` for travel
of body `k` along x and y, per meter, and for its rotation about z through
`centers[k]`, per radian. Every conductor is its own body unless `bodies`
gives conductors that move together one index, and a body rotates about the
centroid of its conductors unless `centers` says otherwise. With potentials
`v`, `0.5 * v @ gradient[:, :, k, 0] @ v` is the force on body `k` along x
per meter of thickness. `hessian[a, b, k, i, l, j]` is the second derivative
for motion `i` of body `k` and motion `j` of body `l`; contracted with
potentials, `-0.5 * v @ hessian[:, :, k, i, l, j] @ v` is the electrostatic
stiffness between those two motions. Both are exact derivatives of the
computed capacitance, so forces on all bodies sum to zero. The point a body
rotates about travels with it.

A pair solved by the field requires `element_size`. Outlines are cut into
straight elements of at most that length, shorter toward convex corners by up to `corner_grading`.
The computed capacitance never exceeds the exact one and rises toward it as
`element_size` shrinks, so halving it shows how far a value has settled. A
mesh of at most `max_elements` is one dense solution; a larger one is solved
in windows, as for the capacitance model, and `asymmetry` then measures what
the windows leave out. Formulas are in
[sidewall capacitance](theory.md#sidewall-capacitance) and
[sidewalls as facing plates](theory.md#sidewalls-as-facing-plates).

## Solve

```python
result = structure.modal()
result = structure.modal(cf.SolveOptions(n_modes=10, method="sparse"))

result.frequencies_hz  # ascending
result.modes  # one column per mode, mass-normalized; rows follow result.dofs
result.dofs  # (body index, "ux" | "uy" | "uz" | "tx" | "ty" | "tz") per row
result.mass, result.stiffness  # assembled 6N by 6N matrices, anchor removed
result.to_npz("modes.npz")

displacement = structure.static({1: [1e-6, 0, 0, 0, 0, 0]})  # fx, fy, fz, tx, ty, tz
```

Matrices are assembled from body mass matrices and effective beam stiffness,
equations (33) to (35), with the anchor's rows and columns dropped. The dense
solver returns all modes or the lowest `n_modes`; the sparse solver uses
shift-invert iteration about `sigma` and needs `n_modes` below the number of
degrees of freedom. Static loads act at the named body's center of mass.

Mode shapes mix translations in meters and rotations in radians, so compare
components of the same kind when reading them.

## Static equilibrium

```python
model = structure.capacitance(params)
drive = structure.electrode_at((1.2e-3, 0.4e-3), "device").net
result = structure.equilibrium(model, {drive: 40.0})

result.displacement  # one value per degree of freedom in result.dofs, from rest
result.stable  # whether a small push is met by a restoring force
result.margin  # share of beam stiffness left where the field softens most: 1 to 0
result.softest  # squared angular frequency of the softest motion under bias
result.stiffness  # tangent stiffness: mechanical plus electrostatic
result.force  # electrostatic forces at the equilibrium
result.potentials, result.charge  # per electrode
result.model  # capacitance model at the equilibrium
result.model.clearance, result.model.closest  # narrowest gap against its width at rest, and where
result.loose  # bodies that no beam holds, kept at rest
result.floating  # nets that held a fixed charge

result = structure.equilibrium(model, {drive: 40.0}, loads={1: [0, 0, -1e-6, 0, 0, 0]})
result = structure.equilibrium(model, {drive: 40.0}, charges={sense: 0.0})  # sense floats
result = structure.equilibrium(model, {drive: 45.0}, start=result)  # from the last equilibrium
```

`structure.equilibrium` finds where the structure rests at one set of
potentials on its nets. It needs both analyses: the structure for beam
stiffness and a capacitance model of that structure for the forces.
Potentials are in volts by net index, as `structure.nets` and
`electrode.net` number them. `charges` names nets that float and the charge
each holds, in coulombs: such a net is connected to nothing, its potential
follows from its charge and moves as bodies move, and `result.potentials`
reports where it ends up. A net in neither is at zero volts, so a net that
floats in the device has to be named in `charges`, with 0.0 if it was never
charged. `loads` adds mechanical forces, as for `static`.

The solve is a Newton iteration that evaluates the capacitance model at
every displacement it visits, so forces and the tangent stiffness are those
of the displaced structure and not a linearization about rest. It starts at
rest, or at `start`: an earlier result or a displacement of every degree of
freedom. A structure under voltage often has more than one equilibrium, and
the iteration finds the one near its start, so the start chooses. `stable`
tells whether that equilibrium holds against a small push, and `margin` is
how far it is from losing that: the share of beam stiffness that the field
has left in the direction it softens most, 1 without voltage, 0 at pull-in,
and negative for an equilibrium that tips over. `softest` is the squared
angular frequency of the softest motion under bias. It also reaches zero at
pull-in, but the motion that gives way is often a stiff one, such as a
plate moving toward the handle, and the softest frequency can stay almost
unchanged until just before. `result.model.clearance` is the narrowest gap
between two bodies as a fraction of its width at rest; nothing in the model
describes bodies that touch, so how small a gap to accept is yours to
decide. A body that no chain of beams joins to the anchor cannot be in
equilibrium under a force; it stays at rest and is listed in `loose`. When
the iteration finds nothing, `cf.NoEquilibriumError` is raised.

The library solves one operating point per call. Sequences are yours to
write, and each solve should start from the equilibrium of the one before
at nearby potentials: that follows the state a device reaches as its
voltages change, where a single solve from rest at a high voltage may land
in another equilibrium. A ramp, a pull-in search, and the phases of a
stepper are a few lines each:

```python
# ramp: displacement against voltage
state, curve = None, []
for volts in np.linspace(5.0, 50.0, 10):
    state = structure.equilibrium(model, {drive: volts}, start=state)
    curve.append((volts, state.displacement, state.margin))

# pull-in: bisect between the last stable equilibrium and the first voltage without one
low, high, state = 0.0, 200.0, None
while high - low > 0.01:
    try:
        trial = structure.equilibrium(model, {drive: (low + high) / 2}, start=state)
    except cf.NoEquilibriumError:
        trial = None
    if trial is not None and trial.stable and trial.model.clearance > 0.1:
        low, state = (low + high) / 2, trial
    else:
        high = (low + high) / 2

# stepper: energize phases in turn, each step from where the last left the rotor
state = None
for phase in [phase_a, phase_b, phase_c] * 3:  # each a list of nets
    state = structure.equilibrium(model, dict.fromkeys(phase, 60.0), start=state)
```

Each of these is quasi-static: it assumes the potentials change slowly
against the motion of the bodies. A drive that switches in a fraction of a
period sets bodies ringing, and where they settle then depends on inertia
and damping. If `plates_cutoff` is lowered so far that teeth standing tip
to tip are taken as plates, the iteration fails at any voltage; that is the
cusp of the plates method and not pull-in. Formulas are in
[static equilibrium](theory.md#static-equilibrium).

## Small signal

```python
from pathlib import Path

bias = structure.equilibrium(model, {shuttle: 40.0})
drive = structure.electrode_at((1.2e-3, 0.4e-3), "device")
sense = structure.electrode_at((0.2e-3, 0.4e-3), "device")
signal = structure.small_signal(bias, [drive, sense], cf.SmallSignalOptions(quality_factor=5000.0))

signal.frequencies_hz  # resonant frequencies under bias, ascending
signal.modes  # mass-normalized mode shapes under bias; rows follow signal.dofs
signal.transduction  # mode by node: force on each mode per volt on each node
signal.capacitance, signal.conductance  # node by node, floating nets eliminated
signal.nodes  # electrodes of each node: ports first, then inner nodes

hertz = np.linspace(50e3, 60e3, 2001)
admittance = signal.admittance(hertz)  # frequencies by ports by ports, complex siemens
motion = signal.displacement(hertz, [0.1, 0.0])  # complex amplitude of every degree of freedom
inductance, capacitance, resistance = signal.motional()  # mode by node
check = signal.linearity([0.1, 0.0], 55e3)  # check.voltage, check.gap, check.peak
Path("resonator.cir").write_text(signal.netlist("resonator"))
```

`structure.small_signal` linearizes the structure about an equilibrium for
small potentials at its ports. Port order is the order of rows and columns
of `admittance` and of drive amplitudes. The equilibrium has to be stable.
The quality factor applies to every mode and is required, since nothing in
the library computes damping.

### Ports and nets

`equilibrium` takes potentials by net; `small_signal` takes ports as
electrodes. A net is everything joined electrically, which shares one
potential when no current flows. An electrode is a piece of a net that
moves with one body, and links that may have resistance join the
electrodes of a net. The two are numbered apart, in `structure.nets` and
`structure.electrodes`. A port is where a source connects: one electrode,
which the source holds at the potential its net has in the equilibrium
plus the signal. Ports are `Electrode` objects, as `electrode_at` returns
them for a point on the pad the source is bonded to; an integer is refused
with a `TypeError`, so a net index cannot pass for a port.

```python
pad = structure.electrode_at((1.2e-3, 0.4e-3), "device")  # the port
pad.net  # its net, the key of its potential in equilibrium
[e for e in structure.electrodes if e.net == pad.net]  # every electrode of that net
```

The equilibrium and the ports give each net one of three roles:

| Net in the equilibrium | In the small-signal model |
| --- | --- |
| Has a potential and a port on one of its electrodes | The source holds the port. Electrodes that links without resistance join to it are the same node; electrodes behind resistance are inner nodes, reached through the resistor network. |
| Has a potential and no port | Held at its potential on every electrode, a signal ground. Its current is not reported. A net named nowhere in the equilibrium is at zero volts and is one of these. |
| Holds a charge | One node that keeps its charge, so its potential follows the signal and the motion. It cannot carry a port, since no source holds it. |

So:

- Without `Material.resistivity`, beams have no resistance, every net is
  one node, and a port stands for its whole net: any electrode of the net
  gives the same model.
- With resistivity, the electrode matters. Choose the pad the source is
  bonded to; a shuttle behind long thin beams then sees its drive late and
  weakened at high frequency.
- A wire with resistance leaves the pad at its far end as an inner node of
  the same net; an ideal wire makes the two pads one node.
- Two ports on one net need resistance between them. Two ports on
  electrodes that links without resistance join would be one node held by
  two sources, which is an error.
- To see the current a net at a fixed potential draws, such as a sense
  electrode at zero volts, give it a port and a drive of zero. Its row of
  `admittance` is the current through that source per volt at each port.
- `signal.nodes` lists the electrode indices of each node, ports first, and
  `structure.electrodes[a].net` the net of each.

### Response

`admittance` gives, at each frequency you pass in hertz, the current into
each port per volt at each port with the others held; it is symmetric.
`displacement` gives complex amplitudes of every degree of freedom for a
drive, one complex amplitude per port, so `abs` is the amplitude and
`np.angle` the phase against the drive. `motional` gives the series branch
of each mode at each node with every other node held: for a single port
near one resonance, `signal.capacitance` in parallel with that branch is the
Butterworth-Van Dyke circuit. Frequencies are yours to choose; for a
resonance, sample a few times more densely than its width
`signal.frequencies_hz[n] / quality_factor`.

### Range of the linearization

The linearization holds while the drive is small against the bias across
each moving gap and the motion small against the gap. `linearity` measures
both for a drive at one frequency: `voltage` is the largest swing across a
pair of electrodes whose capacitance moves, against its bias, and `gap` how
far the narrowest gap closes at the peak of the motion, against its width
in the equilibrium. It warns with `cf.LinearizationWarning` when either
exceeds `linear_tolerance`, and evaluates the capacitance model twice. Use
`warnings.simplefilter("error", cf.LinearizationWarning)` to make that an
error.

### Netlist

`netlist` returns a SPICE subcircuit as text, with pins `port0`, `port1`,
... in port order and node 0 as the potential of nets held without a port:
capacitors and resistors between nodes and one series branch per mode,
coupled to the nodes through controlled sources. Every mode is included
unless `modes` lists some; a mode left out also takes its compliance below
resonance with it. Formulas are in [small signal](theory.md#small-signal).

## Progress and logging

Each stage of `extract`, `modal`, `static`, `plates`, `capacitance`,
`sidewall_capacitance`, `equilibrium`, and `small_signal` logs when it starts and when it finishes, with elapsed time and counts, through the standard `connflex_core`
logger. The logger is silent until a handler is attached:

```python
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
```

```text
2026-11-03 12:00:01,452 heal: started
2026-11-03 12:00:13,850 heal: done in 12.40 s, device 1203, handle 7
```

For a progress bar, pass any object with `start(stage, total)`, `advance(n)`,
and `finish()` as `progress`. `total` is the number of `advance` calls to
expect, or `None` when a stage is one indivisible operation. With tqdm:

```python
from tqdm import tqdm


class Bar:
    def start(self, stage, total):
        self.bar = tqdm(desc=stage, total=total)

    def advance(self, n=1):
        self.bar.update(n)

    def finish(self):
        self.bar.close()


structure = cf.extract(layout, stack, params, progress=Bar())
result = structure.modal(progress=Bar())
```

`extract` reports heal, release, prune, nets, opening, bodies, candidates,
beams, mass, and electrodes, in that order. `modal` and `static` report
assemble and solve, and `equilibrium` reports assemble and equilibrium. Prune and bodies advance once per polygon, beams once per
candidate, mass once per body, and electrodes once per beam, rejected
candidate, and via; the other stages are single operations. `plates` reports
one stage, plates, that advances once per device electrode.
`sidewall_capacitance` with method `"field"` reports mesh, assemble, solve,
derivatives, and curvature; assemble advances once per boundary element and
derivatives once per pair of bodies. For a mesh solved in windows it reports
mesh and windows, which advances once per window. With method `"plates"` it
reports one stage, facing. `capacitance` reports plates, clusters, and
sidewalls; sidewalls advances once per cluster, and once per window for a
cluster solved in windows.
