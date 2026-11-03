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

Polygons within one role are unioned, invalid polygons are rejected with
role and index named, and `layout.device` and friends are shapely
`MultiPolygon` objects in meters.

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

## Process files

Stack, parameters, and solve options load from one TOML file with the same
field names, so a process description can be versioned next to a layout. See
the [parameters page](parameters.md#toml-files) for the format.

```python
stack, params, options = cf.load_config("soi_25um.toml")
```

## Extract the structure

```python
structure = cf.extract(layout, stack, params)

structure.bodies  # rigid bodies, anchor first; polygons per layer in meters
structure.bodies[1].inertial  # mass, center of mass, inertia tensor, 6 by 6 mass matrix
structure.beams  # beams joining exactly two bodies: length, width, thickness, end points
structure.rejected  # beam candidates that did not qualify, each with a reason
structure.graph  # networkx MultiGraph, bodies as nodes and beams as edges
structure.nets  # electrical nets of remaining device and handle polygons
structure.to_dict()  # JSON-serializable tree with GeoJSON-style polygons
```

Stages follow the [theory of operation](theory.md): device and handle layers
are healed and simplified, release holes are filled and used with the
`release` layer to estimate surviving buried oxide, an anchor handle polygon
is selected, polygons with no connection to it through the stack are dropped,
bodies are assembled by overlap through the opened device layer, oxide, and
handle, and device regions removed by the opening become beams when they
touch exactly two bodies and their width lies within `beam_width_min` and
`beam_width_max`. Each body carries its mass, center of mass, and inertia
tensor, computed layer by layer from the stack with exact triangulated
quadrature and with release openings removed from the device layer;
`body.inertial.matrix` is the 6 by 6 mass matrix in body degrees of freedom.

A candidate that touches one body is a stub, one that touches three or more
is a junction; both appear in `structure.rejected` with that reason, which is
the first place to look when an expected beam is missing.

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
