# Input conventions

## Layer roles

| Role | Required | A polygon means |
| --- | --- | --- |
| `device` | yes | device-layer silicon is present here |
| `handle` | yes | handle-layer silicon is present here |
| `release` | no | a release opening is etched here, so device silicon is absent |
| `vias` | no | the device layer is electrically connected to the handle layer here |

## Positive polygons

Every polygon marks where the feature named by its layer is present. No layer
is interpreted as an inverted mask, and no dark-field or clear-field semantics
are applied. The `release` role follows the same rule: a release polygon says
there is a release opening here.

Readers check that polygons are non-empty and valid. Overlapping or touching
polygons within one role are unioned.

## Optional layers

- **`release` absent:** holes inside device polygons with area at or below
  `release_hole_max_area` are the release holes. They enter the oxide
  survival estimate and are then filled. When a `release` layer is given it is
  unioned with such holes.
- **`vias` absent:** no vias. Each connected device island and each connected
  handle island is its own net.

## Electrical interpretation

- Connected `device` polygons form one net; connected `handle` polygons form
  one net.
- A `vias` polygon merges every device net and every handle net it overlaps.
- A via overlapping no conductor on one side produces a warning and is
  ignored.
- Vias are electrical only; the mechanical model ignores them.

## Units

SI throughout: meters, square meters, kilograms, pascals, hertz. The GDS
reader applies the unit stored in the file. The shapely reader takes a
`unit` scale factor.
