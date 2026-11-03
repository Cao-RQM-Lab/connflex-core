# Usage

All quantities are SI. The GDS reader applies the unit stored in the file;
the other readers take a `unit` scale factor, for example `unit=1e-6` for
polygons drawn in micrometers.

## Load a layout

A layout has four polygon roles: `device`, `handle`, and optionally `release`
and `vias`. See the [input conventions](conventions.md).

=== "GDS file"

    ```python
    import connflex_core as cf

    layout = cf.load_gds(
        "resonator.gds",
        layers=cf.LayerMap(device=(1, 0), handle=(2, 0), release=(3, 0), vias=(4, 0)),
    )
    ```

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

=== "gdsfactory"

    ```python
    import connflex_core as cf
    from connflex_core.adapters.gdsfactory_ import from_component

    layout = from_component(component, layers=cf.LayerMap(device=(1, 0), handle=(2, 0)))
    ```

=== "klayout"

    ```python
    import connflex_core as cf
    from connflex_core.adapters.klayout_ import from_cell

    layout = from_cell(cell, layers=cf.LayerMap(device=(1, 0), handle=(2, 0)))
    ```

## Describe the process

```python
stack = cf.Stack(
    device=cf.Layer(thickness=25e-6, material=cf.Material.silicon()),
    box=cf.Layer(thickness=2e-6, material=cf.Material.oxide()),
    handle=cf.Layer(thickness=400e-6, material=cf.Material.silicon()),
)
```

## Choose extraction parameters

Three values have no default because they are process and design specific.
The rest are listed on the [parameters page](parameters.md).

```python
params = cf.ExtractionParams(
    release_hole_max_area=100e-12,
    release_distance=5e-6,
    body_min_width=10e-6,
)
```

Stack, parameters, and solve options can also be loaded from one TOML file:

```python
stack, params, options = cf.load_config("soi_25um.toml")
```

## Extract

```python
structure = cf.extract(layout, stack, params)

structure.bodies  # mass, center of mass, inertia, polygons per layer
structure.beams  # length, width, thickness, endpoints, bodies joined
structure.anchor  # the immovable body claimed from the largest handle polygon
structure.graph  # networkx.Graph, bodies as nodes, beams as edges
structure.to_tree()  # nested dict, JSON serializable
```

## Solve

```python
result = structure.modal()
result = structure.modal(cf.SolveOptions(n_modes=10, method="sparse"))

result.frequencies_hz  # ascending
result.modes  # one column per mode, 6 DOF per body
result.M, result.K  # assembled mass and stiffness matrices
result.dof_index  # maps (body, dof) to a row of M and K
result.to_npz("modes.npz")

u = structure.static({structure.bodies[1]: cf.Load(fz=1e-6)})
```

## Electrical nets

Connected `device` polygons form one net, connected `handle` polygons form
one net, and a `vias` polygon merges the device and handle nets it overlaps.

```python
nets = structure.nets
nets.graph  # conductors as nodes, vias as edges
nets.net_of(beam)
nets.nets_of(body)  # a body can carry several nets
```
