from collections.abc import Callable

import numpy as np
import pytest
from shapely import Polygon, box
from shapely.affinity import rotate

from connflex_core import (
    ElectrostaticParams,
    ExtractionParams,
    Layer,
    Layout,
    Material,
    Stack,
    Structure,
    extract,
)

UM = 1e-6
EPS0 = 8.8541878188e-12
GAP = 2 * UM
ANCHOR = box(0, 0, 400, 400)
BEAM = box(400, 195, 600, 205)
PLATE = box(600, 150, 700, 250)
HANDLE = box(-50, -50, 750, 450)
Params = Callable[..., ExtractionParams]
LayoutUm = Callable[..., Layout]


@pytest.fixture
def shuttle(stack: Stack, params: Params, layout_um: LayoutUm) -> Structure:
    """Bonded pad, beam, and a released 100 um square plate, all over one handle."""
    layout = layout_um([ANCHOR, BEAM, PLATE], [HANDLE])
    return extract(layout, stack, params(release_distance=60 * UM))


@pytest.fixture
def overhang(stack: Stack, params: Params, layout_um: LayoutUm) -> Structure:
    """The same plate half over a handle that ends at x = 650 um, a release hole on that edge."""
    hole = box(645, 220, 655, 228)
    perforated = Polygon(PLATE.exterior, [hole.exterior])
    layout = layout_um([ANCHOR, BEAM, perforated], [box(-50, -50, 650, 450)])
    return extract(layout, stack, params(release_distance=60 * UM))


@pytest.fixture
def turned(stack: Stack, params: Params, layout_um: LayoutUm) -> Structure:
    """A released 500 by 100 um plate on the beam, the whole layout turned by 45 degrees."""
    long_plate = box(600, 150, 1100, 250)
    parts = [rotate(part, 45, origin=(0, 0)) for part in (ANCHOR, BEAM, long_plate)]
    handle = rotate(box(-50, -50, 1150, 450), 45, origin=(0, 0))
    return extract(layout_um(parts, [handle]), stack, params(release_distance=60 * UM))


def test_lowest_gap_is_at_a_corner_of_the_overlap_at_any_orientation(turned: Structure) -> None:
    """Tilted about its own long axis, the plate closes its gap by half its width times the
    tilt at the far long edge, whichever way the layout is drawn; the corners of the
    bounding box of the overlap lie far outside the plate."""
    (plate,) = [p for p in turned.plates() if p.bodies == (1, 0)]
    along = np.array([np.cos(np.pi / 4), np.sin(np.pi / 4)])
    for theta in (2e-2, -2e-2):
        tilt = [0, 0, 0, theta * along[0], theta * along[1], 0]
        assert plate.lowest(device=tilt) == pytest.approx(GAP - 50 * UM * abs(theta), rel=1e-9)
        value, _, _ = plate.evaluate(device=tilt)
        assert value > plate.released
    closing = [0, 0, 0, 4.4e-2 * along[0], 4.4e-2 * along[1], 0]
    assert plate.lowest(device=closing) < 0
    with pytest.raises(ValueError, match="closes"):
        plate.evaluate(device=closing)


def test_released_plate_over_the_handle_is_a_parallel_plate(shuttle: Structure) -> None:
    plate = shuttle.plates()[1]
    assert (plate.index, plate.electrodes, plate.bodies) == (1, (1, 2), (1, 0))
    area = (100 * UM) ** 2
    assert plate.oxide == 0.0
    assert plate.released == pytest.approx(EPS0 * area / GAP)
    value, gradient, hessian = plate.evaluate()
    assert value == pytest.approx(EPS0 * area / GAP)
    assert gradient[2] == pytest.approx(-EPS0 * area / GAP**2)
    assert gradient[8] == pytest.approx(EPS0 * area / GAP**2)
    assert hessian[2, 2] == pytest.approx(2 * EPS0 * area / GAP**3)

    closer = 2 * GAP / 3
    value, gradient, hessian = plate.evaluate(device=[0, 0, -GAP / 3, 0, 0, 0])
    assert value == pytest.approx(EPS0 * area / closer)
    assert gradient[2] == pytest.approx(-EPS0 * area / closer**2)
    assert hessian[2, 2] == pytest.approx(2 * EPS0 * area / closer**3)
    assert plate.evaluate(handle=[0, 0, GAP / 3, 0, 0, 0])[0] == pytest.approx(value)
    with pytest.raises(ValueError, match="gap"):
        plate.evaluate(device=[0, 0, -GAP, 0, 0, 0])

    assert plate.extent == pytest.approx([50 * UM, 50 * UM])
    assert plate.lowest() == GAP
    tilted = plate.lowest(device=[0, 0, -0.5 * UM, 1e-3, 0, 0])
    assert tilted == pytest.approx(GAP - 0.5 * UM - 1e-3 * 50 * UM)


def test_tilt_follows_the_exact_integral_to_second_order(shuttle: Structure) -> None:
    plate = shuttle.plates()[1]
    side = 100 * UM
    _, gradient, hessian = plate.evaluate()
    assert hessian[3, 3] == pytest.approx(EPS0 * side**4 / (6 * GAP**3))
    assert hessian[4, 4] == pytest.approx(EPS0 * side**4 / (6 * GAP**3))
    assert gradient[3:5] == pytest.approx([0, 0], abs=1e-9 * EPS0 * side**3 / GAP**2)

    # Handle turned about y through its own center: gap is GAP + theta (x - cx) under the plate.
    theta = 1e-3
    near, far = (np.array([600, 700]) * UM) - shuttle.anchor.inertial.center[0]
    low, high = GAP + theta * near, GAP + theta * far
    exact = EPS0 * side / theta * np.log(high / low)
    slope = EPS0 * side * (-np.log(high / low) / theta**2 + (far / high - near / low) / theta)
    value, gradient, _ = plate.evaluate(handle=[0, 0, 0, 0, theta, 0])
    assert value == pytest.approx(exact, rel=1e-6)  # second order in tilt is 1.6e-4 of it
    assert gradient[10] == pytest.approx(slope, rel=1e-5)


def test_oxide_under_a_pad_is_a_fixed_capacitor(shuttle: Structure) -> None:
    pad = shuttle.plates()[0]
    assert (pad.electrodes, pad.bodies) == ((0, 2), (0, 0))
    bonded = shuttle.anchor.box.area
    assert bonded == pytest.approx((280 * UM) ** 2, rel=1e-3)
    assert pad.oxide == pytest.approx(3.9 * EPS0 * bonded / GAP)
    assert pad.released == pytest.approx(EPS0 * ((400 * UM) ** 2 - bonded) / GAP)

    moved = np.array([3e-7, -2e-7, 1e-7, 1e-4, -2e-4, 3e-4])
    value, gradient, _ = pad.evaluate(device=moved, handle=moved)
    assert value == pad.oxide + pad.released
    assert np.all(gradient[:6] + gradient[6:] == 0)

    humid = shuttle.plates(ElectrostaticParams(ambient_permittivity=2.0))[0]
    assert humid.oxide == pad.oxide
    assert humid.released == pytest.approx(2 * pad.released)


def test_release_openings_are_removed_from_the_released_plate(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    holes = [box(610, 160, 618, 168), box(640, 200, 648, 208)]
    drawn = box(690, 240, 710, 245)
    perforated = Polygon(PLATE.exterior, [hole.exterior for hole in holes])
    layout = layout_um([ANCHOR, BEAM, perforated], [HANDLE], release=[drawn])
    plate = extract(layout, stack, params(release_distance=60 * UM)).plates()[1]
    area = (100 * 100 - 2 * 8 * 8 - 10 * 5) * UM * UM
    assert plate.released == pytest.approx(EPS0 * area / GAP)


def test_overlap_follows_in_plane_travel_where_the_handle_outline_crosses(
    overhang: Structure,
) -> None:
    plate = overhang.plates()[1]
    assert plate.released == pytest.approx(EPS0 * (50 * 100 - 5 * 8) * UM**2 / GAP)
    cx, cy, cz = overhang.bodies[1].inertial.center
    assert cx == pytest.approx(650 * UM)
    crossing = (100 - 8) * UM  # handle edge under silicon
    lever = ((250**2 - 150**2) / 2 - (228**2 - 220**2) / 2) * UM**2 - crossing * cy
    first = -(50 * 100 * 25 - 5 * 8 * 2.5) * UM**3  # released area times x from the center
    _, gradient, hessian = plate.evaluate()
    scale = EPS0 * crossing / GAP
    assert gradient[0] == pytest.approx(-scale)
    assert gradient[6] == pytest.approx(scale)
    assert gradient[1] == pytest.approx(0, abs=1e-9 * scale)
    assert gradient[5] == pytest.approx(EPS0 * lever / GAP)
    assert gradient[4] == pytest.approx(EPS0 * first / GAP**2 - scale * (-26 * UM - cz))
    assert hessian[0, 2] == pytest.approx(scale / GAP)
    assert hessian[0, 3] == pytest.approx(EPS0 * lever / GAP**2)


def test_gradient_and_hessian_are_derivatives_of_the_capacitance(overhang: Structure) -> None:
    plate = overhang.plates()[1]
    scale = np.array([UM, UM, 0.1 * UM, 1e-3, 1e-3, 1e-3] * 2)
    state = np.random.default_rng(0).uniform(-1, 1, 12) * scale
    value, gradient, hessian = plate.evaluate(state[:6], state[6:])
    step = 1e-4
    for k in range(12):
        shift = np.zeros(12)
        shift[k] = step * scale[k]
        up = plate.evaluate((state + shift)[:6], (state + shift)[6:])
        down = plate.evaluate((state - shift)[:6], (state - shift)[6:])
        slope = (up[0] - down[0]) / (2 * step)
        assert slope == pytest.approx(gradient[k] * scale[k], rel=1e-6, abs=1e-9 * value)
        curvature = (up[1] - down[1]) * scale / (2 * step)
        expected = hessian[k] * scale * scale[k]
        assert curvature == pytest.approx(expected, rel=1e-6, abs=1e-9 * value)


def test_only_electrostatic_analysis_needs_the_dielectric_constant_of_the_box(
    stack: Stack, params: Params, layout_um: LayoutUm
) -> None:
    glass = Layer(thickness=GAP, material=Material.oxide(relative_permittivity=None))
    bare = Stack(device=stack.device, box=glass, handle=stack.handle)
    layout = layout_um([ANCHOR, BEAM, PLATE], [HANDLE])
    structure = extract(layout, bare, params(release_distance=60 * UM))
    assert structure.modal().frequencies_hz.size == 6
    with pytest.raises(ValueError, match="relative_permittivity"):
        structure.plates()
