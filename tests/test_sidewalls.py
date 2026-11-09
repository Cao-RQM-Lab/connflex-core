import warnings
from dataclasses import replace
from typing import Any, cast

import numpy as np
import pytest
from numpy.typing import NDArray
from shapely import MultiPolygon, Polygon, affinity, box, union_all

from connflex_core import ElectrostaticParams, SidewallCapacitance, sidewall_capacitance
from connflex_core.electrostatics import windows
from connflex_core.electrostatics.dense import Charges, charge_derivatives, solve_charges
from connflex_core.electrostatics.outlines import Mesh, mesh_outlines
from connflex_core.progress import quiet

EPS0 = 8.8541878188e-12
UM = 1e-6


def circle(radius: float, sides: int, x: float = 0.0, y: float = 0.0) -> Polygon:
    angle = 2 * np.pi * (np.arange(sides) + 0.5) / sides
    return Polygon(np.column_stack([x + radius * np.cos(angle), y + radius * np.sin(angle)]))


def um(polygon: Polygon) -> Polygon:
    exterior = np.asarray(polygon.exterior.coords) * UM
    return Polygon(exterior, [np.asarray(ring.coords) * UM for ring in polygon.interiors])


def params(size: float, grading: float = 16.0, limit: int = 8000) -> ElectrostaticParams:
    """Field solved whole for every pair, as far as ``limit`` allows."""
    return ElectrostaticParams(
        element_size=size * UM,
        corner_grading=grading,
        max_elements=limit,
        interaction_distance=1e-3,
        dense_cutoff=float("inf"),
        plates_cutoff=float("inf"),
    )


def test_eccentric_circles_have_the_closed_form_capacitance() -> None:
    a, b, offset = 1.0, 2.5, 0.6
    inner = um(circle(a, 180, x=offset))
    outer = um(Polygon(circle(4.0, 180).exterior, [circle(b, 180).exterior]))
    result = sidewall_capacitance([inner, outer], params(1.0))
    expected = 2 * np.pi * EPS0 / np.arccosh((a * a + b * b - offset * offset) / (2 * a * b))
    c = result.capacitance
    assert c[0, 0] == pytest.approx(expected, rel=1e-4)
    assert c == pytest.approx(np.array([[1, -1], [-1, 1]]) * c[0, 0], rel=1e-12)
    outside = np.hypot(*(result.mesh.start / UM).T) > 3.0
    assert outside.sum() == 180
    assert abs(result.charge[outside, 0].sum()) < 1e-5 * c[0, 0]

    ratio = (a * a + b * b - offset * offset) / (2 * a * b)
    slope = expected / np.arccosh(ratio) / np.sqrt(ratio * ratio - 1) * offset / (a * b * UM)
    assert result.gradient[0, 0, 0, 0] == pytest.approx(slope, rel=1e-3)
    assert result.gradient[0, 0, 1, 0] == pytest.approx(-result.gradient[0, 0, 0, 0], rel=1e-9)
    assert result.gradient[0, 0, 0, 1:] == pytest.approx([0, 0], abs=1e-9 * slope)


def test_long_parallel_strips_approach_the_parallel_plate() -> None:
    gap = 2.0

    def strips(length: float) -> list[Polygon]:
        return [um(box(0, 0, length, 4)), um(box(0, 4 + gap, length, 8 + gap))]

    short = sidewall_capacitance(strips(80.0), params(1.0))
    long = sidewall_capacitance(strips(160.0), params(1.0))
    added = long.capacitance[0, 0] - short.capacitance[0, 0]
    assert added == pytest.approx(EPS0 * 80.0 / gap, rel=1e-2)
    assert long.capacitance[0, 0] > EPS0 * 160.0 / gap

    pull = long.gradient[0, 0, 1, 1] - short.gradient[0, 0, 1, 1]
    assert pull == pytest.approx(-EPS0 * 80.0 / (gap * gap * UM), rel=1e-3)

    middle = (long.mesh.start + long.mesh.end) / 2 / UM
    facing = np.argmin(np.hypot(middle[:, 0] - 80.0, middle[:, 1] - 4.0))
    length = np.hypot(*(long.mesh.end[facing] - long.mesh.start[facing]))
    assert long.mesh.conductor[facing] == 0
    assert long.charge[facing, 0] / length == pytest.approx(EPS0 / (gap * UM), rel=1e-6)


def test_refining_the_mesh_converges_from_below() -> None:
    squares = [um(box(0, 0, 10, 10)), um(box(12, 0, 22, 10))]
    coarse, medium, fine = (
        sidewall_capacitance(squares, params(size)).capacitance[0, 0] for size in (2.0, 1.0, 0.5)
    )
    assert coarse < medium < fine
    assert fine - medium < 0.5 * (medium - coarse)
    uniform = sidewall_capacitance(squares, params(1.0, grading=1.0)).capacitance[0, 0]
    assert fine - medium < 0.2 * (fine - uniform)


def test_gradient_and_hessian_are_derivatives_of_the_capacitance() -> None:
    shapes = [box(0, 0, 10, 4), box(12.3, 1, 20, 9), Polygon([(3, 6.1), (9, 7), (8, 12), (2, 11)])]
    bodies = [0, 1, 0]
    centers = np.array([[4.0, 5.0], [15.0, 4.0]])

    def solve(body: int = 0, dof: int = 0, step: float = 0.0) -> SidewallCapacitance:
        moved, pivots = [], centers.copy()
        if dof < 2:
            pivots[body, dof] += step  # a body turns about a point that travels with it
        for shape, owner in zip(shapes, bodies, strict=True):
            if owner == body and dof < 2:
                shape = affinity.translate(shape, *np.eye(2)[dof] * step)
            elif owner == body:
                shape = affinity.rotate(shape, step, origin=tuple(centers[body]), use_radians=True)
            moved.append(um(shape))
        return sidewall_capacitance(moved, params(0.83), bodies=bodies, centers=pivots * UM)

    rest = solve()
    assert rest.gradient.shape == (3, 3, 2, 3)
    assert rest.hessian.shape == (3, 3, 2, 3, 2, 3)
    step = 1e-5
    for body in range(2):
        for dof in range(3):
            ahead, behind = solve(body, dof, step), solve(body, dof, -step)
            per = 2 * step * (UM if dof < 2 else 1.0)
            slope = (ahead.capacitance - behind.capacitance) / per
            expected = rest.gradient[:, :, body, dof]
            assert expected == pytest.approx(slope, rel=1e-6, abs=1e-7 * np.max(abs(expected)))
            bend = (ahead.gradient - behind.gradient) / per
            for other in range(2):
                for travel in range(3):
                    expected = rest.hessian[:, :, other, travel, body, dof]
                    found = bend[:, :, other, travel]
                    assert expected == pytest.approx(
                        found, rel=1e-5, abs=1e-6 * np.max(abs(expected))
                    )


def test_hessian_between_two_cylinders_matches_the_closed_form() -> None:
    distance = 3.0
    result = sidewall_capacitance(
        [um(circle(1.0, 360)), um(circle(1.0, 360, x=distance))], params(1.0)
    )
    spread = np.arccosh(distance / 2)
    root = np.sqrt(distance * distance / 4 - 1)
    slope = -np.pi * EPS0 / (2 * spread**2 * root) / UM
    bend = (
        np.pi
        * EPS0
        / 4
        * (2 / (spread**3 * root**2) + distance / 2 / (spread**2 * root**3))
        / UM**2
    )
    curvature = result.hessian[0, 0]
    assert curvature[1, 0, 1, 0] == pytest.approx(bend, rel=1e-3)
    assert curvature[0, 0, 1, 0] == pytest.approx(-bend, rel=1e-3)
    assert curvature[1, 1, 1, 1] == pytest.approx(slope / (distance * UM), rel=1e-3)
    assert curvature[1, 0, 1, 1] == pytest.approx(0, abs=1e-9 * bend)


def test_moving_every_body_together_changes_nothing() -> None:
    shapes = [um(box(0, 0, 10, 4)), um(box(12.3, 1, 20, 9)), um(box(2, 6, 9, 12))]
    apart = sidewall_capacitance(shapes, params(1.0))
    centroids = np.array([shape.centroid.coords[0] for shape in shapes])
    assert apart.centers == pytest.approx(centroids)
    pivot = np.tile([[1.0 * UM, 2.0 * UM]], (3, 1))
    together = sidewall_capacitance(shapes, params(1.0), centers=pivot)
    gradient, curvature = together.gradient, together.hessian
    scale = abs(gradient).max(axis=(0, 1, 2))
    assert gradient == pytest.approx(gradient.transpose(1, 0, 2, 3), abs=1e-12 * np.max(scale))
    assert np.all(abs(gradient.sum(axis=2)).max(axis=(0, 1)) < 1e-9 * scale)
    assert np.all(abs(gradient.sum(axis=0)).max(axis=(0, 1)) < 1e-9 * scale)
    bend = np.max(abs(curvature[..., :2, :, :2]))
    assert curvature == pytest.approx(curvature.transpose(0, 1, 4, 5, 2, 3), abs=1e-9 * bend)
    assert np.max(abs(curvature[..., :2].sum(axis=4)[..., :2, :])) < 1e-8 * bend


def interleaved(pitch: float = 14.0, fingers: int = 6) -> list[Polygon]:
    """Two blocks whose fingers interleave with gaps of 3 um."""
    left = [box(0, pitch * k, 60, pitch * k + 4) for k in range(fingers + 1)]
    right = [box(20, pitch * k + 7, 80, pitch * k + 11) for k in range(fingers)]
    top = pitch * fingers + 4
    return [
        um(cast(Polygon, union_all([box(-30, 0, 0, top), *left]))),
        um(cast(Polygon, union_all([box(80, 0, 110, top), *right]))),
    ]


def in_windows(conductors: list[Polygon]) -> ElectrostaticParams:
    """Windows of an eighth of the mesh, which may hold half of it."""
    regions = [MultiPolygon([conductor]) for conductor in conductors]
    count = len(mesh_outlines(regions, 3 * UM, 4.0).conductor)
    return replace(
        params(3.0, grading=4.0, limit=count // 2),
        window_elements=count // 8,
        interaction_distance=30 * UM,
    )


def test_mesh_above_the_element_limit_is_solved_in_windows() -> None:
    """Force along the fingers is what windows find least well: it comes from finger
    tips, which face open space."""
    conductors = interleaved()
    whole = sidewall_capacitance(conductors, params(3.0, grading=4.0))
    windowed = sidewall_capacitance(conductors, in_windows(conductors))
    assert whole.windows == 1
    assert windowed.windows >= 8
    assert (whole.method[0, 1], windowed.method[0, 1]) == ("dense", "windowed")
    assert whole.asymmetry < 1e-9 < windowed.asymmetry < 2e-2
    for found, expected in (
        (windowed.capacitance, whole.capacitance),
        (windowed.gradient, whole.gradient),
        (windowed.hessian, whole.hessian),
    ):
        assert found == pytest.approx(expected, abs=3e-2 * np.max(abs(expected)))
        assert np.array_equal(found, found.swapaxes(0, 1))
        assert np.max(abs(found.sum(axis=1))) < 1e-12 * np.max(abs(found))
    assert windowed.charge.sum(axis=0) == pytest.approx(0.0, abs=1e-2 * windowed.capacitance[0, 0])


def test_window_that_is_a_translate_takes_the_solution_of_its_original(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A comb and its copy 400 um off are cut alike, so the copy's windows take the
    solutions of the comb's. Its bodies turn about points that are not translates, and
    its derivatives still match the copy solved alone."""
    comb = interleaved()
    copy = [affinity.translate(conductor, 400 * UM) for conductor in comb]
    limited = in_windows(comb)
    solved = []

    def counted(*args: Any, **kwargs: Any) -> Charges:
        solved.append(1)
        return solve_charges(*args, **kwargs)

    monkeypatch.setattr(windows, "solve_charges", counted)
    pivots = np.array([[-15.0, 40.0], [95.0, 50.0]]) * UM
    alone = sidewall_capacitance(copy, limited, centers=pivots)
    once = len(solved)
    both = sidewall_capacitance([*comb, *copy], limited, centers=np.vstack([pivots, pivots]))
    assert len(solved) == 2 * once
    assert both.windows == 2 * alone.windows
    copied = np.arange(2, 4)
    for found, expected in (
        (both.capacitance[np.ix_(copied, copied)], alone.capacitance),
        (both.gradient[np.ix_(copied, copied, copied)], alone.gradient),
        (both.hessian[np.ix_(copied, copied, copied, range(3), copied)], alone.hessian),
    ):
        assert found == pytest.approx(expected, rel=1e-9, abs=1e-12 * np.max(abs(expected)))


def test_field_of_a_window_is_solved_once_however_its_search_ends(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A window whose overlap search settles before the widest overlap keeps the
    derivatives of the solve that settled it, instead of solving the same field once
    more for them."""
    conductors = interleaved()
    original = windows.solve_window
    solved: list[bytes] = []

    def counted(mesh: Mesh, *rest: Any) -> Any:
        solved.append(mesh.start.tobytes() + rest[1].tobytes())
        return original(mesh, *rest)

    monkeypatch.setattr(windows, "solve_window", counted)
    found = sidewall_capacitance(conductors, in_windows(conductors), workers=1)
    assert found.windows >= 8
    assert len(solved) == len(set(solved))


def test_derivatives_of_charge_on_some_elements_are_exact() -> None:
    """A window counts charge on its core alone; elements keep their places as bodies move."""
    squares = [MultiPolygon([um(box(0, 0, 10, 10))]), MultiPolygon([um(box(12, -2, 20, 9))])]
    mesh = mesh_outlines(squares, 2 * UM, 4.0)
    core = (mesh.start[:, 0] + mesh.end[:, 0]) / 2 > 6 * UM
    body = np.array([0, 1])
    centers = np.array([[5.0, 5.0], [16.0, 3.0]]) * UM

    def solve(dof: int = 0, step: float = 0.0) -> tuple[NDArray[np.float64], ...]:
        pivots, moving = centers.copy(), mesh.conductor == 1
        ends = [mesh.start.copy(), mesh.end.copy()]
        for points in ends:
            if dof < 2:
                points[moving, dof] += step
            else:
                arm = points[moving] - centers[1]
                turn = np.array([[np.cos(step), np.sin(step)], [-np.sin(step), np.cos(step)]])
                points[moving] = centers[1] + arm @ turn
        if dof < 2:
            pivots[1, dof] += step
        charges = solve_charges(Mesh(ends[0], ends[1], mesh.conductor), 2, quiet, core)
        return charges.capacitance, *charge_derivatives(charges, body, pivots, quiet)

    value, gradient, hessian = solve()
    assert not np.allclose(value, value.T, rtol=1e-3)
    for dof, step in ((0, 1e-4 * UM), (1, 1e-4 * UM), (2, 1e-5)):
        ahead, behind = solve(dof, step), solve(dof, -step)
        slope = (ahead[0] - behind[0]) / (2 * step)
        assert gradient[:, :, 1, dof] == pytest.approx(
            slope, rel=1e-6, abs=1e-7 * np.max(abs(slope))
        )
        bend = (ahead[1] - behind[1]) / (2 * step)
        for other in range(2):
            for travel in range(3):
                expected = hessian[:, :, other, travel, 1, dof]
                found = bend[:, :, other, travel]
                assert expected == pytest.approx(found, rel=1e-5, abs=1e-6 * np.max(abs(found)))


def test_each_pair_takes_the_method_its_feature_size_falls_in() -> None:
    """A block faces a long wall across 1.5 um, feature size 133, and a toothed rack
    whose teeth are 2 gaps wide. The rack moves."""
    block = box(0, 0, 200, 30)
    wall = box(0, 31.5, 200, 51.5)
    rack = union_all(
        [box(0, -40, 200, -10)] + [box(20 + 40 * k, -10, 24 + 40 * k, -2) for k in range(5)]
    )
    conductors = [um(block), um(wall), um(cast(Polygon, rack))]
    mixed = ElectrostaticParams(
        element_size=2 * UM, corner_grading=4.0, interaction_distance=30 * UM
    )
    assert (mixed.dense_cutoff, mixed.plates_cutoff) == (10.0, 100.0)
    bodies, centers = [0, 1, 2], np.array([[100, 15], [100, 41.5], [100, -25]]) * UM
    result = sidewall_capacitance(conductors, mixed, bodies=bodies, centers=centers)
    assert result.method[0, 1] == "plates"
    assert result.method[0, 2] == "dense"
    assert result.method[1, 2] == ""
    assert result.feature_size[0, 1] == pytest.approx(200 / 1.5, rel=1e-6)
    assert result.feature_size[0, 2] < 10
    c = result.capacitance
    assert np.array_equal(c, c.T)
    assert np.max(abs(c.sum(axis=1))) < 1e-12 * np.max(abs(c))
    plates = sidewall_capacitance(conductors, replace(mixed, dense_cutoff=0.0, plates_cutoff=0.0))
    dense = sidewall_capacitance(
        conductors, replace(mixed, dense_cutoff=float("inf"), plates_cutoff=float("inf"))
    )
    assert c[0, 1] == pytest.approx(plates.capacitance[0, 1], rel=1e-12)
    assert c[0, 2] == pytest.approx(dense.capacitance[0, 2], rel=2e-2)
    assert result.gradient[0, 2, 2] == pytest.approx(
        dense.gradient[0, 2, 2], rel=2e-2, abs=2e-2 * abs(dense.gradient[0, 2, 2]).max()
    )
    assert not c[1, 2]

    windowed = sidewall_capacitance(conductors, replace(mixed, plates_cutoff=float("inf")))
    assert windowed.method[0, 1] == "windowed"
    assert windowed.capacitance[0, 1] == pytest.approx(dense.capacitance[0, 1], rel=2e-2)

    def moved(dof: int, step: float) -> SidewallCapacitance:
        shifted = conductors[:2]
        pivots = centers.copy()
        if dof < 2:
            pivots[2, dof] += step
            shifted.append(affinity.translate(conductors[2], *np.eye(2)[dof] * step))
        else:
            origin = tuple(centers[2])
            shifted.append(affinity.rotate(conductors[2], step, origin=origin, use_radians=True))
        return sidewall_capacitance(shifted, mixed, bodies=bodies, centers=pivots)

    for dof, step in ((0, 1e-10), (1, 1e-10), (2, 1e-7)):
        ahead, behind = moved(dof, step), moved(dof, -step)
        slope = (ahead.capacitance - behind.capacitance) / (2 * step)
        expected = result.gradient[:, :, 2, dof]
        assert expected == pytest.approx(slope, rel=1e-4, abs=1e-5 * np.max(abs(expected)))
        bend = (ahead.gradient[:, :, 2, dof] - behind.gradient[:, :, 2, dof]) / (2 * step)
        expected = result.hessian[:, :, 2, dof, 2, dof]
        assert expected == pytest.approx(bend, rel=1e-3, abs=1e-4 * np.max(abs(expected)))


def test_cavity_that_holds_no_conductor_is_not_meshed() -> None:
    solid, other = um(box(0, 0, 10, 10)), um(box(12, 0, 22, 10))
    hollow = Polygon(solid.exterior, [um(box(3, 3, 7, 7)).exterior])
    filled = sidewall_capacitance([solid, other], params(1.0))
    perforated = sidewall_capacitance([hollow, other], params(1.0))
    assert len(perforated.mesh.conductor) == len(filled.mesh.conductor)
    assert np.array_equal(perforated.capacitance, filled.capacitance)


def test_elements_tile_a_long_outline_within_the_size_limit() -> None:
    strip = MultiPolygon([um(box(0, 0, 3000, 7))])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mesh = mesh_outlines([strip], 2.0 * UM, 16.0)
    lengths = np.hypot(*(mesh.end - mesh.start).T)
    assert lengths.sum() == pytest.approx(strip.length)
    assert lengths.max() <= 2.0 * UM
    assert lengths.min() == pytest.approx(2.0 * UM / 16 * np.expm1(1.0), rel=0.2)


def test_inputs_that_cannot_be_solved_are_rejected() -> None:
    squares = [um(box(0, 0, 10, 10)), um(box(12, 0, 22, 10))]
    with pytest.raises(ValueError, match="element_size"):
        sidewall_capacitance(squares, replace(params(1.0), element_size=None))
    with pytest.raises(ValueError, match="max_elements"):
        sidewall_capacitance(squares, params(1.0, limit=50))
    with pytest.raises(ValueError, match="conductors 0 and 1"):
        sidewall_capacitance([squares[0], um(box(10, 0, 20, 10))], params(1.0))
    with pytest.raises(ValueError, match="bodies"):
        sidewall_capacitance(squares, params(1.0), bodies=[0, 2])
    with pytest.raises(ValueError, match="centers"):
        sidewall_capacitance(squares, params(1.0), centers=[[0.0, 0.0]])


def test_vertices_a_rounding_apart_make_no_element() -> None:
    """Healing can leave two vertices a rounding apart; an element that short is noise."""
    sliver = Polygon([(0, 0), (10, 0), (10 + 1e-12, 1e-12), (10, 10), (0, 10)])
    mesh = mesh_outlines([MultiPolygon([um(sliver)])], 2 * UM, 4.0)
    assert np.hypot(*(mesh.end - mesh.start).T).min() > 1e-3 * UM
    square = mesh_outlines([MultiPolygon([um(box(0, 0, 10, 10))])], 2 * UM, 4.0)
    assert len(mesh.conductor) == len(square.conductor)


def test_assembly_evaluates_each_pair_of_elements_once(monkeypatch: pytest.MonkeyPatch) -> None:
    """The matrix of equation (64) is symmetric, so assembly evaluates the pair integral
    of each unordered pair once and mirrors it."""
    from connflex_core.electrostatics import dense
    from connflex_core.electrostatics.pairs import pair_integrals as original

    squares = [MultiPolygon([box(0, 0, 1e-5, 1e-5)]), MultiPolygon([box(1.2e-5, 0, 2.2e-5, 1e-5)])]
    mesh = mesh_outlines(squares, 2.5e-6, 1.0)
    count = len(mesh.conductor)
    sizes: list[int] = []

    def counted(*ends: NDArray[np.complex128]) -> NDArray[np.float64]:
        result = original(*ends)
        sizes.append(result.size)
        return result

    monkeypatch.setattr(dense, "pair_integrals", counted)
    charges = solve_charges(mesh, 2, quiet)
    assert count > 8
    assert sum(sizes) == count * (count + 1) // 2
    assert charges.capacitance[0, 1] == pytest.approx(charges.capacitance[1, 0], rel=1e-12)
