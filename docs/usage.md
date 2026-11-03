# Usage

All quantities are SI. The GDS reader applies the unit stored in the file;
the shapely reader takes a `unit` scale factor, for example `unit=1e-6` for
polygons drawn in micrometers.

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

    The file's unit is applied. Pass `cell="NAME"` when the file has more than
    one top-level cell. Holes stored as cut lines are repaired into polygons
    with holes.

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

Polygons within one role are unioned, invalid polygons are rejected with the
role and index named, and `layout.device` and friends are shapely
`MultiPolygon` objects in meters.

Readers are looked up by name through `cf.adapters.get(name)`; `cf.adapters.available()`
lists them. A package adds a reader with `cf.adapters.register(name, function)` or by
publishing an entry point in the `connflex_core.adapters` group.

## Describe the process

The stack lists the three SOI layers from the top down. Thickness is required
for each; the material presets can be overridden field by field.

```python
stack = cf.Stack(
    device=cf.Layer(thickness=25e-6, material=cf.Material.silicon()),
    box=cf.Layer(thickness=2e-6, material=cf.Material.oxide()),
    handle=cf.Layer(thickness=400e-6, material=cf.Material.silicon(youngs_modulus=130e9)),
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
field names, so a process description can be versioned next to a layout. The
file format is on the [parameters page](parameters.md#toml-files).

```python
stack, params, options = cf.load_config("soi_25um.toml")
```
