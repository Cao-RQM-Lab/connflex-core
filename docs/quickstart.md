# Quick start

Three short guides, each a complete script on one small design: a square plate
with release holes, hung on four beams between two pads that stay bonded to
the handle. The first takes the layout to its resonant modes, the second finds
its capacitance to the handle and the force that a voltage puts on it, and the
third solves the equilibrium under that voltage, linearizes about it, and
follows a ramp in time. Each guide continues from the one before, so run them
in order. The [usage page](usage.md) covers every option in depth.

Everything is in SI units. Polygons below are drawn in micrometers and
converted with `unit=1e-6`; a GDS file is read with `cf.load_gds` and needs no
unit, as the [usage page](usage.md#load-a-layout) shows.

## 1. From a layout to resonant modes

A layout has two required roles, `device` and `handle`, each a list of
polygons that mark where that silicon is present. The plate carries a grid of
small holes, so that the release etch undercuts it completely while the wide
pads keep oxide under their middles and stay anchored.

```python
import numpy as np
from shapely import Polygon, box

import connflex_core as cf

holes = [box(c - 2, r - 2, c + 2, r + 2) for c in range(10, 100, 20) for r in range(10, 100, 20)]
plate = Polygon(box(0, 0, 100, 100).exterior, [hole.exterior for hole in holes])
beams = [box(x, y, x + 200, y + 10) for x in (-200, 100) for y in (10, 80)]
pads = [box(-500, -100, -200, 200), box(300, -100, 600, 200)]
layout = cf.Layout.from_shapely(
    device=[plate, *beams, *pads], handle=[box(-600, -200, 700, 300)], unit=1e-6
)
```

The process is a stack of three layers with a thickness each, and three
extraction parameters have no default because they depend on the process and
the design: the largest hole that counts as a release opening, the distance
the release etch undercuts, and the narrowest feature that is a rigid body
rather than a beam.

```python
stack = cf.Stack(
    device=cf.Layer(thickness=25e-6, material=cf.Material.silicon()),
    box=cf.Layer(thickness=2e-6, material=cf.Material.oxide()),
    handle=cf.Layer(thickness=400e-6, material=cf.Material.silicon()),
)
params = cf.ExtractionParams(
    release_hole_max_area=100e-12,  # holes of 4 um are release openings
    release_distance=12e-6,  # the etch undercuts 12 um from every opening and edge
    body_min_width=20e-6,  # anything narrower is a beam candidate
)
structure = cf.extract(layout, stack, params)
```

The structure holds the rigid bodies, the anchor first, and the beams between
them. Everything that was not accepted as a beam is listed with its reason,
which is the first place to look when a flexure is missing.

```python
len(structure.bodies), len(structure.beams)  # (2, 4): the anchor with both pads, the plate, four beams
structure.bodies[1].inertial.mass  # kilograms, release holes removed
structure.beams[0].length, structure.beams[0].width  # meters
structure.rejected  # candidates that are not beams, each with a reason
```

Modal analysis assembles the mass and stiffness matrices of the bodies and
beams and solves for the resonances. The result lists frequencies in hertz,
the mass-normalized mode shapes, and the degree of freedom of each row.

```python
modes = structure.modal()
modes.frequencies_hz  # ascending; six modes for one movable body
modes.dofs[2]  # (1, "uz"): the plate's travel toward the handle
up_and_down = int(np.argmax(np.abs(modes.modes[2])))  # the mode that moves it most that way
f_vertical = modes.frequencies_hz[up_and_down]
```

## 2. Capacitance and electrostatic force

The plate over the handle is a capacitor across the released gap, and the
pads are capacitors through their oxide. The capacitance model collects every
such plate, and with sidewalls off needs no further parameter. Rows of the
matrix sum to zero, so the capacitance between two electrodes is the negative
of their off-diagonal entry.

```python
model = structure.capacitance(cf.ElectrostaticParams(include_sidewalls=False))
plate_electrode = structure.electrode_at((50e-6, 50e-6), "device")
handle_electrode = structure.electrode_at((50e-6, 50e-6), "handle")
gap_capacitance = -model.matrix[plate_electrode.index, handle_electrode.index]  # farads
```

With potentials on the electrodes the model gives the generalized force on
every degree of freedom and the electrostatic stiffness, both exact derivatives
of the capacitance. Potentials are per electrode; one value per net spreads
over its electrodes through `model.net`.

```python
handle_net = handle_electrode.net
volts = np.where(model.net == handle_net, 10.0, 0.0)  # 10 V on the handle, pads and plate at 0 V
force = model.force(volts)
uz = model.dofs.index((plate_electrode.body, "uz"))
force[uz]  # negative: the plate is pulled toward the handle
softening = model.stiffness(volts)[uz, uz]  # negative: the field softens the vertical spring
```

Capacitance between sidewalls in the layout plane comes from a field solution
on the outlines. The standalone function takes polygons in meters and returns
capacitance per meter of thickness, with its derivatives for the motion of
each conductor.

```python
wall_a = box(0.0, 0.0, 100e-6, 10e-6)
wall_b = box(0.0, 13e-6, 100e-6, 23e-6)  # a 3 um trench between two walls 100 um long
walls = cf.sidewall_capacitance(
    [wall_a, wall_b], cf.ElectrostaticParams(element_size=1e-6, interaction_distance=10e-6)
)
-walls.capacitance[0, 1] * 25e-6  # farads for a 25 um layer, a little above eps * 100 / 3 * 25 um
walls.method[0, 1]  # how the pair was solved, chosen by the length of the walls over their gap
```

## 3. Equilibrium, small signal, and transient

Under a voltage the plate moves until its beams balance the pull. The solve
takes potentials by net, evaluates the capacitance model at every
displacement it visits, and reports the equilibrium with its stability and
the margin left before pull-in.

```python
bias = structure.equilibrium(model, {handle_net: 20.0})
bias.displacement[uz]  # meters toward the handle, negative
bias.stable, bias.margin  # True, and the share of beam stiffness the field has left
bias.model.clearance  # the narrowest gap against its width at rest
```

The small-signal model linearizes about that equilibrium for signals at
ports. A port is an electrode, the pad a source is bonded to; the handle,
which has a potential and no port, is the signal ground. The admittance is
evaluated at the frequencies you choose, and a SPICE netlist of the model is
text to write wherever you like.

```python
pad = structure.electrode_at((-350e-6, 50e-6), "device")
signal = structure.small_signal(bias, [pad], cf.SmallSignalOptions(quality_factor=100.0))
coupled = int(np.argmax(np.abs(signal.transduction[:, 0])))  # the mode the port drives most
f_biased = signal.frequencies_hz[coupled]  # below f_vertical: the field softens the spring
hertz = np.linspace(0.9, 1.1, 401) * f_biased
admittance = signal.admittance(hertz)[:, 0, 0]  # complex siemens into the pad
netlist = signal.netlist("plate")  # SPICE subcircuit as text
```

A transient follows the motion in time under sources. A source drives one
electrode with a waveform, here a ramp to 20 V through the driver's
resistance; the handle is held at zero volts because nothing drives it. The
run ends at the last time given, or earlier if a gap closes.

```python
period = 1 / f_biased
ramp = cf.PiecewiseLinear([0.0, 2 * period], [0.0, 20.0])
source = cf.Source(pad, ramp, resistance=1e3)
times = np.linspace(0.0, 25 * period, 501)
run = structure.transient(model, [source], times, cf.TransientOptions(quality_factor=20.0))
run.displacement[-1, uz]  # settles at the equilibrium of the same voltage
run.current[-1, 0]  # amperes the source drives at the end, near zero at rest
run.energy.imbalance[-1]  # what the integration left of the energy balance, in joules
```

## Where to go next

- The [usage page](usage.md) walks through every input and result, with
  sequences such as voltage ramps and pull-in searches.
- The [parameters page](parameters.md) lists every field with its default and
  its effect.
- The [theory pages](theory/index.md) give the equations behind each number,
  and the [assumptions page](assumptions.md) says what each idealization
  costs.
