import numpy as np
import pytest
from shapely import Point, Polygon, affinity, box

from connflex_core import ElectrostaticParams, SidewallCapacitance, sidewall_capacitance

EPS0 = 8.8541878188e-12
UM = 1e-6
PLATES = ElectrostaticParams(sidewall_method="plates", interaction_distance=30 * UM)


def um(polygon: Polygon) -> Polygon:
    return Polygon(np.asarray(polygon.exterior.coords) * UM)


def test_offset_strips_are_one_plate_over_their_overlap() -> None:
    gap, overlap = 2.0, 90.0
    lower, upper = um(box(0, 0, 100, 4)), um(box(10, 4 + gap, 110, 8 + gap))
    result = sidewall_capacitance([lower, upper], PLATES)
    c = result.capacitance
    assert c[0, 1] == pytest.approx(-EPS0 * overlap / gap, rel=1e-12)
    assert c == pytest.approx(np.array([[1, -1], [-1, 1]]) * c[0, 0], rel=1e-12)
    assert result.gradient[0, 0, 1, :2] == pytest.approx(
        [-EPS0 / (gap * UM), -EPS0 * overlap / (gap**2 * UM)], rel=1e-9
    )
    assert result.gradient[0, 0, 0, :2] == pytest.approx(-result.gradient[0, 0, 1, :2], rel=1e-9)
    assert result.hessian[0, 0, 1, 1, 1, 1] == pytest.approx(
        2 * EPS0 * overlap / (gap**3 * UM**2), rel=1e-9
    )
    assert result.hessian[0, 0, 1, 0, 1, 1] == pytest.approx(EPS0 / (gap**2 * UM**2), rel=1e-9)
    assert len(result.charge) == len(result.mesh.conductor) == 0
    walls = result.walls
    assert sorted(zip(walls.conductor.tolist(), walls.facing.tolist(), strict=True)) == [
        (0, 1),
        (1, 0),
    ]
    assert np.sort(np.concatenate([walls.start[:, 0], walls.end[:, 0]])) == pytest.approx(
        np.array([10, 10, 100, 100]) * UM
    )
    assert walls.gap == pytest.approx(gap * UM)


def test_a_wall_sees_only_the_nearest_conductor_within_reach() -> None:
    stack = [um(box(0, 6 * k, 50, 6 * k + 4)) for k in range(3)]
    c = sidewall_capacitance(stack, PLATES).capacitance
    assert c[0, 1] == pytest.approx(-EPS0 * 50 / 2)
    assert c[1, 2] == pytest.approx(-EPS0 * 50 / 2)
    assert c[0, 2] == 0.0
    far = [um(box(0, 0, 50, 4)), um(box(0, 40, 50, 44))]
    assert not sidewall_capacitance(far, PLATES).capacitance.any()


def test_wall_faced_lies_to_the_right_of_each_piece() -> None:
    """A block in the opening of a frame: walls of a hole face inward."""
    hole = np.array([(0, 0), (50, 0), (50, 14), (0, 14)])
    frame = Polygon(np.array([(-10, -10), (60, -10), (60, 24), (-10, 24)]) * UM, [hole * UM])
    conductors = [frame, um(box(3, 2, 47, 11))]
    walls = sidewall_capacitance(conductors, PLATES).walls
    assert set(zip(walls.conductor.tolist(), walls.facing.tolist(), strict=True)) == {
        (0, 1),
        (1, 0),
    }
    along = walls.end - walls.start
    right = np.stack([along[:, 1], -along[:, 0]], axis=1) / np.hypot(*along.T)[:, None]
    for point, faced in zip(walls.start + right * walls.gap[:, :1], walls.facing, strict=True):
        assert conductors[faced].boundary.distance(Point(point)) < 1e-9 * UM


def test_tilted_wall_integrates_the_inverse_gap() -> None:
    lower = um(box(0, 0, 40, 4))
    upper = um(Polygon([(0, 6), (40, 10), (40, 14), (0, 14)]))
    c = sidewall_capacitance([lower, upper], PLATES).capacitance
    from_below = 10 * np.log(6 / 2)  # gap 2 + x / 10 over 0 to 40
    last = (40 - 0.2) / 1.01  # where the normal of the tilted wall still lands on the strip
    from_above = 10 * np.log((2 + last / 10) / 2)
    assert -c[0, 1] == pytest.approx(EPS0 * (from_below + from_above) / 2, rel=1e-9)


def test_flank_that_leans_by_a_rounding_error_is_no_plate() -> None:
    """Two teeth tip to tip. One flank is square to the tips or leans across the corner
    of the other tooth by 0.01 nm, as coordinates of a turned layout do."""
    tooth = um(box(40, 0, 60, 4))
    square, leaning = (
        sidewall_capacitance(
            [tooth, um(Polygon([(40, 6), (60 - lean, 6), (60 + lean, 11), (40, 11)]))], PLATES
        )
        for lean in (0.0, 5e-6)
    )
    assert -square.capacitance[0, 1] == pytest.approx(EPS0 * 20 / 2, rel=1e-12)
    assert leaning.capacitance == pytest.approx(square.capacitance, rel=1e-6)
    plate = 2 * EPS0 * 20 / (2**3 * UM**2)
    assert leaning.hessian[0, 1, 1, 1, 1, 1] == pytest.approx(-plate, rel=1e-6)
    assert abs(leaning.hessian[0, 1, 1, 0, 1, 0]) < 1e-6 * plate


def test_gradient_and_hessian_are_derivatives_of_the_plate_capacitance() -> None:
    """A block between two walls hides part of each from the other, so its motion
    changes capacitance between them too. No wall lies along the normal of another:
    there capacitance has a kink in rotation."""
    shapes = [
        box(0, 0, 30, 4),
        Polygon([(-2, 10), (33, 10.6), (33, 14), (-2, 14)]),
        Polygon([(10, 6), (18, 6.3), (17.5, 8), (10.4, 8.2)]),
    ]
    centers = np.array([[4.0, 5.0], [15.0, 12.0], [12.0, 3.0]])

    def solve(body: int = 0, dof: int = 0, step: float = 0.0) -> SidewallCapacitance:
        moved, pivots = list(shapes), centers.copy()
        if dof < 2:
            pivots[body, dof] += step
            moved[body] = affinity.translate(shapes[body], *np.eye(2)[dof] * step)
        else:
            origin = tuple(centers[body])
            moved[body] = affinity.rotate(shapes[body], step, origin=origin, use_radians=True)
        return sidewall_capacitance([um(shape) for shape in moved], PLATES, centers=pivots * UM)

    rest = solve()
    assert np.all(rest.capacitance[np.triu_indices(3, 1)] < 0)
    assert np.any(rest.gradient[0, 1, 2])
    step = 1e-5
    for body in range(3):
        for dof in range(3):
            ahead, behind = solve(body, dof, step), solve(body, dof, -step)
            per = 2 * step * (UM if dof < 2 else 1.0)
            slope = (ahead.capacitance - behind.capacitance) / per
            expected = rest.gradient[:, :, body, dof]
            assert expected == pytest.approx(slope, rel=1e-6, abs=1e-7 * np.max(abs(expected)))
            bend = (ahead.gradient - behind.gradient) / per
            for other in range(3):
                for travel in range(3):
                    expected = rest.hessian[:, :, other, travel, body, dof]
                    found = bend[:, :, other, travel]
                    assert expected == pytest.approx(
                        found, rel=1e-5, abs=1e-6 * np.max(abs(rest.hessian))
                    )


def test_method_and_its_parameters_must_be_chosen() -> None:
    strips = [um(box(0, 0, 100, 4)), um(box(0, 6, 100, 10))]
    with pytest.raises(ValueError, match="sidewall_method"):
        sidewall_capacitance(strips, ElectrostaticParams(element_size=1 * UM))
    with pytest.raises(ValueError, match="interaction_distance"):
        sidewall_capacitance(strips, ElectrostaticParams(sidewall_method="plates"))
    with pytest.raises(ValueError, match="sidewall_method"):
        ElectrostaticParams(sidewall_method="guess")  # type: ignore[arg-type]
