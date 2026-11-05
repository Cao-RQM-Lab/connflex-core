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

structure.bodies  # rigid bodies, anchor first; polygons per layer in meters
structure.bodies[1].inertial  # mass, center of mass, inertia tensor, 6 by 6 mass matrix
structure.bodies[1].stubs  # device silicon attached to this body alone; mass, no stiffness
structure.beams  # beams joining exactly two bodies: length, width, thickness, end points
structure.rejected  # beam candidates that did not qualify, each with a reason
structure.graph  # networkx MultiGraph, bodies as nodes and beams as edges
structure.nets  # electrical nets of remaining device and handle polygons
structure.electrodes  # silicon with one potential and one rigid motion: nodes of the electrical model
structure.links  # beams, bridges, and vias that conduct between electrodes
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
directly and needs no structure.

## Electrical nets

```python
nets = structure.nets
nets.count  # number of nets
nets.conductors  # connected device regions, then handle regions, each with its net
nets.vias  # each via with the conductors it joins on both sides
nets.graph  # networkx MultiGraph, conductors as nodes and vias as edges
nets.conductor_at((x, y), "device")  # conductor covering a point

structure.beams[0].net  # net of the beam's device conductor
structure.bodies[0].nets  # net of each device polygon of the body
```

Each connected device region and each connected handle region is one
conductor. A via joins every conductor it overlaps on both sides, and nets
are the connected components. A via that overlaps no conductor on one side
is ignored with a warning. Polygons dropped by anchor pruning are not
conductors. `cf.find_nets` computes the same from any lists of polygons.

## Electrodes and links

```python
structure.electrodes  # device electrodes, then handle electrodes
structure.links  # beams, bridges, and vias that conduct between electrodes
structure.electrode_at((x, y), "device")  # electrode covering a point

electrode = structure.electrodes[0]
electrode.body, electrode.conductor, electrode.net
electrode.region  # body polygon followed by its stubs, in meters

link = structure.links[0]
link.electrodes  # indices of the two electrodes joined
link.kind  # "beam", "bridge", or "via"
link.source  # index in structure.beams, structure.rejected, or structure.nets.vias
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
A via links the device and handle electrodes it overlaps. Electrodes joined
through links always share a net.

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

model.matrix  # n by n over structure.electrodes, farads
model.oxide, model.released, model.sidewall  # the three parts it sums
-model.matrix[a, b]  # capacitance between electrodes a and b
model.net_matrix  # between nets, electrodes of each net tied together
model.clusters  # groups of device electrodes whose sidewalls were solved together
model.plates  # plates behind the oxide and released parts

model.dofs  # (body, "ux" | "uy" | "uz" | "tx" | "ty" | "tz") per degree of freedom
model.gradient  # n by n by len(dofs): change of the matrix with each
force = model.force(potentials)  # one potential per electrode
softening = model.stiffness(potentials)  # adds to the mechanical stiffness
```

`structure.capacitance` collects every capacitance between electrodes at
rest into one matrix: plates between device and handle electrodes, through
oxide and across released gaps, and sidewalls between device electrodes. It
needs no mechanical solution. Rows sum to zero, so only differences of
potential matter and nothing outside the structure holds charge.

Sidewalls are solved in clusters. Device electrodes whose silicon lies
within `interaction_distance` of each other, directly or through a chain of
neighbors, form a cluster and are solved together with `element_size`;
electrodes in different clusters share no sidewall capacitance. A cluster
that needs more boundary elements than `max_elements` raises an error. With
`skip_large_clusters=True` it is left out instead and appears in
`model.clusters` with `solved` false, so a model can be built for the rest.

The model also holds how the matrix changes when bodies move, in the six
degrees of freedom of every movable body and in the order of `result.dofs`
from `modal`. `force` returns the generalized forces for given electrode
potentials and `stiffness` the electrostatic stiffness, whose negative
diagonal entries soften the suspension. Potentials are per electrode;
`np.asarray(per_net)[model.net]` spreads one value per net. Sidewall forces
act at mid-height of the device layer, so they also tilt a body whose center
of mass lies lower.

Either part can be left out, and then needs none of its parameters:
`cf.ElectrostaticParams(include_sidewalls=False)` gives plates alone with no
element size, and `include_plates=False` gives sidewalls alone. Definitions
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
params = cf.ElectrostaticParams(element_size=1e-6)
result = cf.sidewall_capacitance([finger, slot], params)

result.capacitance  # n by n, farads per meter of thickness
-result.capacitance[0, 1] * 25e-6  # farads between the two across a 25 um layer
result.gradient  # n by n by bodies by 3: change with x, y, and rotation of each body
result.hessian  # n by n by bodies by 3 by bodies by 3: second derivatives
result.centers  # point each body rotates about
result.charge  # charge on each element per volt on each conductor
result.mesh.start, result.mesh.end, result.mesh.conductor  # boundary elements

result = cf.sidewall_capacitance(conductors, params, bodies=[0, 0, 1], centers=centers)
```

`sidewall_capacitance` solves the field between conductors drawn in the
layout plane. Each conductor is a shapely polygon or multipolygon in meters
and is one equipotential; conductors must not overlap or touch. The field is
two-dimensional, which holds where the layer is thick against its trenches,
so the result is per meter of thickness.

The conductors carry no net charge together and potential far away floats.
Only differences of potential matter: every row of `capacitance` sums to
zero, and `-capacitance[a, b]` is the capacitance between conductors `a` and
`b`.

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

`element_size` is required. Outlines are cut into straight elements of at
most that length, shorter toward convex corners by up to `corner_grading`.
The computed capacitance never exceeds the exact one and rises toward it as
`element_size` shrinks, so halving it shows how far a value has settled. The
solve is dense, and a mesh above `max_elements` raises an error; pass fewer
conductors or a larger `element_size`. Formulas are in
[sidewall capacitance](theory.md#sidewall-capacitance).

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

## Progress and logging

Each stage of `extract`, `modal`, `static`, `plates`, `capacitance`, and
`sidewall_capacitance` logs when it starts and when it finishes, with elapsed time and counts, through the standard `connflex_core`
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
assemble and solve. Prune and bodies advance once per polygon, beams once per
candidate, mass once per body, and electrodes once per beam, rejected
candidate, and via; the other stages are single operations. `plates` reports
one stage, plates, that advances once per device electrode.
`sidewall_capacitance` reports mesh, assemble, solve, derivatives, and
curvature; assemble advances once per boundary element and derivatives once
per pair of bodies. `capacitance` reports plates, clusters, and sidewalls;
sidewalls advances once per cluster.
