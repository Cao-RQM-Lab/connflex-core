# Input conventions

## Layer roles

| Role | Required | A polygon means |
| --- | --- | --- |
| `device` | yes | device-layer silicon is present here |
| `handle` | yes | handle-layer silicon is present here |
| `release` | no | a release opening is etched here, so device silicon is absent |
| `vias` | no | device and handle layers are electrically connected here |

## Positive polygons

Every polygon marks where the feature named by its layer is present. No layer
is interpreted as an inverted mask, and no dark-field or clear-field semantics
are applied. The `release` role follows the same rule: a release polygon says
there is a release opening here.

Readers reject invalid polygons, drop empty ones, and require at least one
polygon on each role given. Overlapping or touching polygons within one role
are unioned.

## Optional layers

- **`release` absent:** holes inside device polygons with area at or below
  `release_hole_max_area` are release holes. They enter the oxide survival
  estimate, are removed from device-layer mass, and are filled for
  connectivity and beam extraction. When a `release` layer is given it is
  unioned with such holes.
- **`vias` absent:** no vias. Each connected device island and each connected
  handle island is its own net.

Release openings, whether drawn on the `release` layer or found as holes,
are modeled as perforations. They enter the oxide estimate and are removed
from device-layer mass, but they never split a body, narrow a beam, or break
a net: a release polygon that severs device silicon leaves one body and one
conductor. The `release` layer has no area threshold and is simplified but
not healed.

## Electrical interpretation

- Connected `device` polygons form one net; connected `handle` polygons form
  one net.
- A `vias` polygon merges every device net and every handle net it overlaps.
- A via overlapping no conductor on one side produces a warning and is
  ignored.
- Each via conducts between the device and handle electrodes it joins with
  `via_resistance`, which is 0 $\Omega$ unless set and does not depend on via area.
- Vias are electrical only; the mechanical model ignores them.

## Units

SI throughout: meters, square meters, kilograms, pascals, hertz. GDS,
KLayout, and gdsfactory readers apply the unit stored in the file or
database. The shapely reader takes a `unit` scale factor.
