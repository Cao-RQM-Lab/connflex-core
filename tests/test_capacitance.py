import json
from dataclasses import replace
from typing import Any

import numpy as np
import pytest
from shapely import MultiPolygon, Polygon, affinity, box

from connflex_core import (
    CapacitanceModel,
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
THICKNESS, BOX = 25 * UM, 2 * UM
GAP, WIDTH, FINGERS, OVERLAP = 3.0, 4.0, 4, 40.0
PITCH = 2 * (GAP + WIDTH)
GAPS = 2 * FINGERS


@pytest.fixture(scope="module")
def comb() -> Structure:
    """Bonded pad with five fingers, and a released shuttle with four between them on two beams."""
    stator = [box(300, 100 + PITCH * k, 360, 100 + PITCH * k + WIDTH) for k in range(FINGERS + 1)]
    rotor = [
        box(320, 100 + PITCH * k + WIDTH + GAP, 380, 100 + PITCH * k + 2 * WIDTH + GAP)
        for k in range(FINGERS)
    ]
    spine, beams = box(380, 60, 420, 200), [box(420, 70, 620, 80), box(420, 180, 620, 190)]
    device = [box(0, 0, 300, 300), *stator, *rotor, spine, *beams, box(620, 0, 920, 300)]
    layout = Layout.from_shapely(device=device, handle=[box(-50, -50, 970, 350)], unit=UM)
    stack = Stack(
        device=Layer(thickness=THICKNESS, material=Material.silicon()),
        box=Layer(thickness=BOX, material=Material.oxide()),
        handle=Layer(thickness=400 * UM, material=Material.silicon()),
    )
    params = ExtractionParams(
        release_hole_max_area=100 * UM * UM, release_distance=60 * UM, body_min_width=20 * UM
    )
    return extract(layout, stack, params)


@pytest.fixture(scope="module")
def settings() -> ElectrostaticParams:
    """Field solved whole for every pair, in the layout plane alone."""
    return ElectrostaticParams(
        element_size=4 * UM,
        interaction_distance=30 * UM,
        dense_cutoff=float("inf"),
        plates_cutoff=float("inf"),
        include_fringe=False,
    )


@pytest.fixture(scope="module")
def model(comb: Structure, settings: ElectrostaticParams) -> CapacitanceModel:
    return comb.capacitance(settings)


@pytest.fixture(scope="module")
def at(comb: Structure) -> dict[str, int]:
    """Electrode indices by name."""
    return {
        "stator": comb.electrode_at((150 * UM, 150 * UM), "device").index,
        "shuttle": comb.electrode_at((400 * UM, 130 * UM), "device").index,
        "pad": comb.electrode_at((770 * UM, 150 * UM), "device").index,
        "handle": comb.electrode_at((0.0, 0.0), "handle").index,
    }


def test_matrix_sums_sidewalls_and_plates_between_electrodes(
    comb: Structure, model: CapacitanceModel, at: dict[str, int]
) -> None:
    stator, shuttle, pad, handle = (at[name] for name in ("stator", "shuttle", "pad", "handle"))
    assert model.matrix == pytest.approx(model.oxide + model.released + model.sidewall)
    assert model.matrix == pytest.approx(model.matrix.T)
    assert np.max(abs(model.matrix.sum(axis=1))) < 1e-9 * np.max(abs(model.matrix))

    ideal = GAPS * EPS0 * THICKNESS * OVERLAP / GAP
    assert ideal < -model.sidewall[stator, shuttle] < 1.3 * ideal
    assert model.sidewall[stator, stator] == pytest.approx(-model.sidewall[stator, shuttle])
    assert not model.sidewall[[pad, handle]].any()
    assert [(cluster.electrodes, cluster.windows) for cluster in model.clusters] == [
        ((stator, shuttle), 1)
    ]

    silicon = (40 * 140 + FINGERS * 60 * WIDTH) * UM**2
    assert model.released[shuttle, handle] == pytest.approx(-EPS0 * silicon / BOX)
    assert model.oxide[shuttle, handle] == 0.0
    assert model.oxide[pad, handle] == pytest.approx(-3.9 * EPS0 * (180 * UM) ** 2 / BOX, rel=1e-2)
    assert model.oxide[handle, handle] == pytest.approx(
        -model.oxide[pad, handle] - model.oxide[stator, handle]
    )


def test_net_matrix_ties_electrodes_of_one_net(
    comb: Structure, model: CapacitanceModel, at: dict[str, int]
) -> None:
    net = {name: comb.electrodes[index].net for name, index in at.items()}
    assert net["pad"] == net["shuttle"] != net["stator"]
    assert model.net.tolist() == [electrode.net for electrode in comb.electrodes]
    nets = model.net_matrix
    assert nets.shape == (comb.nets.count, comb.nets.count)
    assert nets[net["stator"], net["shuttle"]] == pytest.approx(
        model.matrix[at["stator"], at["shuttle"]]
    )
    assert nets[net["shuttle"], net["handle"]] == pytest.approx(
        model.matrix[at["shuttle"], at["handle"]] + model.matrix[at["pad"], at["handle"]]
    )
    assert np.max(abs(nets.sum(axis=1))) < 1e-9 * np.max(abs(nets))


def test_each_part_can_be_left_out_and_needs_only_its_own_parameters(
    comb: Structure, settings: ElectrostaticParams, model: CapacitanceModel
) -> None:
    plates_only = comb.capacitance(ElectrostaticParams(include_sidewalls=False))
    assert not plates_only.sidewall.any()
    assert plates_only.clusters == ()
    assert plates_only.oxide == pytest.approx(model.oxide)

    walls_only = comb.capacitance(replace(settings, include_plates=False))
    assert not walls_only.oxide.any()
    assert not walls_only.released.any()
    assert walls_only.sidewall == pytest.approx(model.sidewall)

    with pytest.raises(ValueError, match="interaction_distance"):
        comb.capacitance(replace(settings, interaction_distance=None))
    with pytest.raises(ValueError, match="element_size"):
        comb.capacitance(replace(settings, element_size=None))

    assert comb.capacitance(replace(settings, interaction_distance=1 * UM)).clusters == ()


def test_default_cutoffs_solve_a_comb_in_windows_and_report_it(
    comb: Structure, model: CapacitanceModel, at: dict[str, int]
) -> None:
    """Fingers overlap 40 um across 3 um gaps: a feature size of about 13 gaps."""
    defaults = comb.capacitance(
        ElectrostaticParams(element_size=4 * UM, interaction_distance=30 * UM)
    )
    (cluster,) = defaults.clusters
    stator, shuttle = (comb.electrodes[at[name]].net for name in ("stator", "shuttle"))
    (interaction,) = cluster.interactions
    assert interaction.nets == tuple(sorted((stator, shuttle)))
    assert interaction.method == "windowed"
    assert 10 < interaction.feature_size < 15
    assert defaults.sidewall == pytest.approx(
        model.sidewall, rel=2e-2, abs=2e-2 * np.max(abs(model.sidewall))
    )
    report = defaults.to_dict()["clusters"][0]["interactions"]
    assert report == [
        {
            "nets": list(interaction.nets),
            "feature_size": interaction.feature_size,
            "method": "windowed",
        }
    ]
    there = defaults.at([1e-6, 0.0, 0.0, 0.0, 0.0, 0.0])
    assert there.clusters[0].interactions == cluster.interactions


def test_cluster_above_the_element_limit_is_solved_in_windows(
    comb: Structure, settings: ElectrostaticParams, model: CapacitanceModel, at: dict[str, int]
) -> None:
    (whole,) = model.clusters
    limited = replace(
        settings, max_elements=whole.elements // 2, window_elements=whole.elements // 8
    )
    windowed = comb.capacitance(limited)
    (cluster,) = windowed.clusters
    assert (cluster.electrodes, cluster.elements) == (whole.electrodes, whole.elements)
    assert cluster.windows >= 8
    assert whole.asymmetry < 1e-9 < cluster.asymmetry < 2e-2
    for found, expected in (
        (windowed.sidewall, model.sidewall),
        (windowed.gradient, model.gradient),
    ):
        assert found == pytest.approx(expected, abs=2e-2 * np.max(abs(expected)))
    potentials = np.zeros(len(comb.electrodes))
    potentials[at["stator"]] = 10.0
    softening = model.stiffness(potentials)
    assert windowed.stiffness(potentials) == pytest.approx(
        softening, abs=2e-2 * np.max(abs(softening))
    )
    place = [1e-6, 0.2e-6, 0.0, 0.0, 0.0, 0.0]
    there, reference = windowed.at(place), model.at(place)
    assert there.clusters[0].windows == cluster.windows
    assert there.sidewall == pytest.approx(
        reference.sidewall, abs=2e-2 * np.max(abs(reference.sidewall))
    )


def test_comb_drives_with_the_pull_of_each_gap_and_the_plate_pulls_down(
    comb: Structure, model: CapacitanceModel, at: dict[str, int]
) -> None:
    assert model.dofs == tuple((1, name) for name in ("ux", "uy", "uz", "tx", "ty", "tz"))
    potentials = np.zeros(len(comb.electrodes))
    potentials[at["stator"]], potentials[at["handle"]] = 10.0, 5.0
    force = model.force(potentials)
    drive = GAPS * 0.5 * EPS0 * THICKNESS * 10.0**2 / (GAP * UM)
    assert force[0] == pytest.approx(-drive, rel=0.03)
    assert abs(force[1]) < 1e-3 * drive
    silicon = (40 * 140 + FINGERS * 60 * WIDTH) * UM**2
    assert force[2] == pytest.approx(-0.5 * EPS0 * silicon * 5.0**2 / BOX**2)

    stiffness = model.stiffness(potentials)
    assert stiffness == pytest.approx(stiffness.T)
    side = GAPS * EPS0 * THICKNESS * OVERLAP * UM * 10.0**2 / (GAP * UM) ** 3
    assert stiffness[1, 1] == pytest.approx(-side, rel=0.06)
    assert stiffness[2, 2] == pytest.approx(-EPS0 * silicon * 5.0**2 / BOX**3)
    with pytest.raises(ValueError, match="potentials"):
        model.force(potentials[:-1])


def test_plates_method_takes_every_comb_gap_as_a_parallel_plate(
    comb: Structure, at: dict[str, int]
) -> None:
    """Finger tips are 20 um from the wall they point at, beyond the 10 um that counts."""
    params = ElectrostaticParams(
        dense_cutoff=0.0,
        plates_cutoff=0.0,
        interaction_distance=10 * UM,
        include_plates=False,
        include_fringe=False,
    )
    model = comb.capacitance(params)
    stator, shuttle = at["stator"], at["shuttle"]
    assert -model.sidewall[stator, shuttle] == pytest.approx(
        GAPS * EPS0 * THICKNESS * OVERLAP / GAP, rel=1e-9
    )
    (cluster,) = model.clusters
    assert (cluster.electrodes, cluster.elements, cluster.windows) == ((stator, shuttle), 0, 0)
    assert len(cluster.walls.conductor) == 16
    assert {interaction.method for interaction in cluster.interactions} == {"plates"}
    assert cluster.walls.gap == pytest.approx(GAP * UM)

    potentials = np.zeros(len(comb.electrodes))
    potentials[stator] = 10.0
    force, stiffness = model.force(potentials), model.stiffness(potentials)
    drive = GAPS * 0.5 * EPS0 * THICKNESS * 10.0**2 / (GAP * UM)
    assert force[:2] == pytest.approx([-drive, 0.0], rel=1e-9, abs=1e-9 * drive)
    side = GAPS * EPS0 * THICKNESS * OVERLAP * UM * 10.0**2 / (GAP * UM) ** 3
    assert stiffness[1, 1] == pytest.approx(-side, rel=1e-9)


def test_model_at_displaced_bodies_follows_the_plate_formulas(
    comb: Structure, at: dict[str, int]
) -> None:
    """Shuttle moved along and across its fingers and toward the handle."""
    rest = comb.capacitance(
        ElectrostaticParams(
            dense_cutoff=0.0, plates_cutoff=0.0, interaction_distance=10 * UM, include_fringe=False
        )
    )
    stator, shuttle, handle = at["stator"], at["shuttle"], at["handle"]
    assert not rest.displacement.any()
    travel, side, rise = 5.0, 1.0, -0.5
    moved = rest.at([travel * UM, side * UM, rise * UM, 0.0, 0.0, 0.0])
    assert moved.displacement.tolist() == [travel * UM, side * UM, rise * UM, 0.0, 0.0, 0.0]
    facing = GAPS / 2 * EPS0 * THICKNESS * (1 / (GAP - side) + 1 / (GAP + side))
    assert -moved.sidewall[stator, shuttle] == pytest.approx(facing * (OVERLAP - travel), rel=1e-9)
    assert moved.gradient[stator, stator, 0] == pytest.approx(-facing / UM, rel=1e-9)
    silicon = (40 * 140 + FINGERS * 60 * WIDTH) * UM**2
    assert -moved.released[shuttle, handle] == pytest.approx(EPS0 * silicon / (BOX + rise * UM))
    assert moved.oxide == pytest.approx(rest.oxide)
    assert rest.clearance == 1.0
    assert moved.clearance == pytest.approx((GAP - side) / GAP)
    assert moved.closest == (stator, shuttle)
    lower = rest.at([0.0, 0.0, -1.0 * UM, 0.0, 0.0, 0.0])
    assert lower.clearance == pytest.approx((BOX - 1.0 * UM) / BOX)
    assert lower.closest == (shuttle, handle)
    with pytest.raises(ValueError, match="displacement"):
        rest.at([0.0, 0.0])


def test_derivatives_at_displaced_bodies_are_slopes_of_the_displaced_model(
    comb: Structure, at: dict[str, int]
) -> None:
    rest = comb.capacitance(
        ElectrostaticParams(
            dense_cutoff=0.0, plates_cutoff=0.0, interaction_distance=30 * UM, include_fringe=False
        )
    )
    place = np.array([2e-6, 0.4e-6, -0.3e-6, 2e-4, -1e-4, 3e-4])
    here = rest.at(place)
    potentials = np.zeros(len(comb.electrodes))
    potentials[at["stator"]], potentials[at["handle"]] = 10.0, 5.0
    stiffness = here.stiffness(potentials)
    for dof in range(6):
        step = np.zeros(6)
        step[dof] = 1e-10 if dof < 3 else 1e-8
        ahead, behind = rest.at(place + step), rest.at(place - step)
        slope = (ahead.matrix - behind.matrix) / (2 * step[dof])
        expected = here.gradient[:, :, dof]
        assert expected == pytest.approx(slope, rel=1e-5, abs=1e-6 * np.max(abs(expected)))
        softening = -(ahead.force(potentials) - behind.force(potentials)) / (2 * step[dof])
        scale = np.sqrt(abs(np.diag(stiffness) * stiffness[dof, dof]))
        assert np.all(abs(stiffness[:, dof] - softening) < 1e-4 * abs(softening) + 1e-5 * scale)


def test_charge_gradient_is_the_slope_of_the_charge(
    model: CapacitanceModel, at: dict[str, int]
) -> None:
    """The charge on each electrode per degree of freedom, from which the force follows."""
    volts = np.zeros(len(model.net))
    volts[at["stator"]] = 10.0
    sensitivity = model.charge_gradient(volts)
    assert sensitivity.shape == (len(model.net), len(model.dofs))
    assert model.force(volts) == pytest.approx(0.5 * volts @ sensitivity)
    for k in (0, 2):  # travel along the fingers, and toward the handle
        step = np.zeros(len(model.dofs))
        step[k] = 1e-10
        slope = (model.at(step).matrix - model.at(-step).matrix) @ volts / (2 * step[k])
        assert sensitivity[:, k] == pytest.approx(slope, rel=1e-6, abs=1e-6 * np.max(abs(slope)))


def test_gradient_and_stiffness_follow_a_moved_structure(
    comb: Structure, settings: ElectrostaticParams, model: CapacitanceModel, at: dict[str, int]
) -> None:
    center = comb.bodies[1].inertial.center
    potentials = np.zeros(len(comb.electrodes))
    potentials[at["stator"]], potentials[at["handle"]] = 10.0, 5.0

    def moved(dof: int, step: float) -> CapacitanceModel:
        shift = np.zeros(3)
        electrodes = []
        for electrode in comb.electrodes:
            region = electrode.region
            if electrode.body == 1 and electrode.layer == "device" and dof < 2:
                shift[dof] = step
                region = affinity.translate(region, *shift[:2])
            elif electrode.body == 1 and electrode.layer == "device":
                region = affinity.rotate(region, step, origin=tuple(center[:2]), use_radians=True)
            electrodes.append(replace(electrode, region=region))
        shuttle = replace(
            comb.bodies[1], inertial=replace(comb.bodies[1].inertial, center=center + shift)
        )
        there = replace(comb, electrodes=tuple(electrodes), bodies=(comb.bodies[0], shuttle))
        return there.capacitance(settings)

    for dof, step in ((0, 1e-6), (5, 1e-3)):
        place = np.zeros(6)
        place[dof] = step
        rebuilt = moved(min(dof, 2), step)
        assert model.at(place).matrix == pytest.approx(
            rebuilt.matrix, rel=1e-7, abs=1e-9 * np.max(abs(model.matrix))
        )

    stiffness = model.stiffness(potentials)
    for dof, step in ((0, 2e-9), (1, 2e-9), (5, 2e-6)):
        ahead, behind = moved(min(dof, 2), step), moved(min(dof, 2), -step)
        slope = (ahead.matrix - behind.matrix) / (2 * step)
        expected = model.gradient[:, :, dof]
        assert expected == pytest.approx(slope, rel=1e-5, abs=1e-6 * np.max(abs(expected)))
        softening = -(ahead.force(potentials) - behind.force(potentials)) / (2 * step)
        for other in (0, 1, 5):
            scale = np.sqrt(abs(stiffness[other, other] * stiffness[dof, dof]))
            assert stiffness[other, dof] == pytest.approx(
                softening[other], rel=1e-4, abs=1e-5 * scale
            )


def toothed_rack() -> Structure:
    """Bar with teeth at 10 um pitch on two long beams, under a fixed block with teeth,
    4 um from tip to tip; the bar reaches 150 um past the block at either end."""
    bar = [box(-300, 0, 300, 30), *(box(x, 30, x + 5, 40) for x in range(-295, 300, 10))]
    block = [box(-150, 54, 150, 200), *(box(x, 44, x + 5, 54) for x in range(-145, 150, 10))]
    beams = [box(-299, -800, -297, 0), box(297, -800, 299, 0)]
    anchors = [box(-450, -1100, -150, -800), box(150, -1100, 450, -800)]
    layout = Layout.from_shapely(
        device=[*bar, *block, *beams, *anchors], handle=[box(-500, -1200, 500, 300)], unit=UM
    )
    stack = Stack(
        device=Layer(thickness=THICKNESS, material=Material.silicon()),
        box=Layer(thickness=BOX, material=Material.oxide()),
        handle=Layer(thickness=400 * UM, material=Material.silicon()),
    )
    params = ExtractionParams(
        release_hole_max_area=100 * UM * UM, release_distance=60 * UM, body_min_width=20 * UM
    )
    return extract(layout, stack, params)


def test_windows_follow_bodies_that_travel_past_their_overlap() -> None:
    """Bar slid up to seven teeth along. Windows kept from rest would miss teeth that now
    face their cores, and the row of teeth would count bar teeth past the end of the
    block as interior: pull 40 percent off. Field beyond ``interaction_distance``, which
    windows and rows drop, leaves a few percent as the overhang of the bar changes."""
    rack = toothed_rack()
    whole_field = ElectrostaticParams(
        element_size=2 * UM,
        interaction_distance=30 * UM,
        dense_cutoff=float("inf"),
        plates_cutoff=float("inf"),
        include_plates=False,
        include_fringe=False,
    )
    whole = rack.capacitance(whole_field)
    (cluster,) = whole.clusters
    windowed = rack.capacitance(
        replace(
            whole_field,
            max_elements=cluster.elements - 1,
            window_elements=cluster.elements // 4,
        )
    )
    assert windowed.clusters[0].windows >= 8
    along = whole.dofs.index((whole.dofs[0][0], "ux"))
    for travel in (12.5, 42.5, 72.5, -32.5):
        place = np.zeros(len(whole.dofs))
        place[along] = travel * UM
        there, reference = windowed.at(place), whole.at(place)
        assert there.sidewall == pytest.approx(
            reference.sidewall, abs=1e-2 * np.max(abs(reference.sidewall))
        )
        pull = reference.gradient[:, :, along]
        assert there.gradient[:, :, along] == pytest.approx(pull, abs=5e-2 * np.max(abs(pull)))


def test_parallel_plate_pulls_and_tilts_a_body_that_carries_handle_silicon() -> None:
    """A fixed pad 4 um from a 200 um plate that is bonded to its own handle island."""
    pad, plate, far = box(0, 0, 300, 300), box(304, 50, 504, 250), box(704, 0, 1004, 300)
    beams = [box(504, 60, 704, 70), box(504, 230, 704, 240)]
    frame = Polygon(box(-50, -50, 1054, 350).exterior, [box(310, 30, 498, 270).exterior])
    layout = Layout.from_shapely(
        device=[pad, plate, *beams, far], handle=[frame, box(330, 70, 478, 230)], unit=UM
    )
    stack = Stack(
        device=Layer(thickness=THICKNESS, material=Material.silicon()),
        box=Layer(thickness=BOX, material=Material.oxide()),
        handle=Layer(thickness=400 * UM, material=Material.silicon()),
    )
    params = ExtractionParams(
        release_hole_max_area=100 * UM * UM, release_distance=60 * UM, body_min_width=20 * UM
    )
    structure = extract(layout, stack, params)
    fixed = structure.electrode_at((150 * UM, 150 * UM), "device").index
    lever = -THICKNESS / 2 - structure.bodies[1].inertial.center[2]
    assert lever > 100 * UM

    model = structure.capacitance(
        ElectrostaticParams(
            element_size=4 * UM,
            interaction_distance=30 * UM,
            dense_cutoff=float("inf"),
            plates_cutoff=float("inf"),
            include_fringe=False,
        )
    )
    potentials = np.zeros(len(structure.electrodes))
    potentials[fixed] = 10.0
    force, stiffness = model.force(potentials), model.stiffness(potentials)
    facing, gap = 200 * UM, 4 * UM
    assert force[0] == pytest.approx(-0.5 * EPS0 * THICKNESS * facing * 10.0**2 / gap**2, rel=0.04)
    assert stiffness[0, 0] == pytest.approx(-EPS0 * THICKNESS * facing * 10.0**2 / gap**3, rel=0.03)
    assert force[4] == pytest.approx(lever * force[0])
    assert stiffness[4, 4] == pytest.approx(lever**2 * stiffness[0, 0])


def test_report_lists_capacitance_by_pair_part_net_and_cluster(
    comb: Structure, model: CapacitanceModel, at: dict[str, int]
) -> None:
    report = json.loads(json.dumps(model.to_dict()))
    pairs = {tuple(pair["electrodes"]): pair for pair in report["pairs"]}
    drive = pairs[tuple(sorted((at["stator"], at["shuttle"])))]
    assert drive["sidewall"] == pytest.approx(-model.sidewall[at["stator"], at["shuttle"]])
    assert drive["fringe"] == 0.0
    assert drive["total"] == drive["sidewall"]
    assert drive["nets"] == [comb.electrodes[i].net for i in drive["electrodes"]]
    (moving,) = drive["gradient"]
    assert moving["body"] == 1
    assert moving["values"] == pytest.approx(
        (-model.gradient[at["stator"], at["shuttle"]]).tolist()
    )
    bonded = pairs[tuple(sorted((at["pad"], at["handle"])))]
    assert bonded["oxide"] > 0
    assert bonded["released"] > 0
    assert bonded["gradient"] == []
    assert np.array(report["net_matrix"]) == pytest.approx(model.net_matrix)
    (cluster,) = model.clusters
    assert report["clusters"] == [
        {
            "electrodes": list(cluster.electrodes),
            "elements": cluster.elements,
            "windows": 1,
            "asymmetry": cluster.asymmetry,
            "anchors": 0,
            "sections": 0,
            "interactions": [
                {
                    "nets": list(interaction.nets),
                    "feature_size": interaction.feature_size,
                    "method": "dense",
                }
                for interaction in cluster.interactions
            ],
        }
    ]
    assert report["unassigned_area"] == pytest.approx(2 * 200 * 10 * UM**2)


def test_workers_share_windows_and_sections_and_change_nothing(comb: Structure) -> None:
    params = ElectrostaticParams(
        element_size=4 * UM, interaction_distance=30 * UM, window_elements=40
    )
    alone = comb.capacitance(params, workers=1)
    shared = comb.capacitance(params, workers=4)
    (cluster,) = shared.clusters
    assert cluster.windows > 1
    assert cluster.sections > 1

    def same(found: np.ndarray, expected: np.ndarray) -> None:
        assert found == pytest.approx(expected, rel=1e-10, abs=1e-12 * np.max(abs(expected)))

    same(shared.matrix, alone.matrix)
    same(shared.gradient, alone.gradient)
    for mine, theirs in zip(shared.blocks, alone.blocks, strict=True):
        same(mine.values, theirs.values)
    place = [1e-6, 0.2e-6, 0.0, 0.0, 0.0, 0.0]
    same(shared.at(place, workers=4).matrix, alone.at(place).matrix)
    with pytest.raises(ValueError, match="workers"):
        comb.capacitance(params, workers=0)


def test_derivatives_about_the_core_corner_are_those_of_a_translate() -> None:
    """A window may take the derivatives of a translate solved in its own coordinates:
    about the corner of each copy's core they are the same, which lets the copy that
    holds the factorization serve whichever copy asks first."""
    from connflex_core.electrostatics.dense import charge_derivatives, solve_charges
    from connflex_core.electrostatics.outlines import mesh_outlines
    from connflex_core.progress import quiet

    def copy(shift: float) -> tuple[np.ndarray, np.ndarray]:
        squares = [box(0, 0, 10, 10), box(13, -1, 21, 8)]
        regions = [affinity.translate(square, shift, 2 * shift) for square in squares]
        scaled = [MultiPolygon([affinity.scale(r, UM, UM, origin=(0, 0))]) for r in regions]
        mesh = mesh_outlines(scaled, 2 * UM, 4.0)
        core = (mesh.start[:, 0] + mesh.end[:, 0]) / 2 > (6 + shift) * UM
        charges = solve_charges(mesh, 2, quiet, core)
        corner = ((mesh.start + mesh.end) / 2)[core].min(axis=0)
        return charge_derivatives(charges, np.array([0, 1]), np.tile(corner, (2, 1)), quiet)

    for mine, theirs in zip(copy(0.0), copy(37.0), strict=True):
        assert mine == pytest.approx(theirs, rel=1e-9, abs=1e-12 * np.max(abs(theirs)))


def test_released_electrode_alone_over_the_handle_keeps_its_fringe() -> None:
    """A device electrode with no neighbor within reach forms no cluster, but with the
    handle beneath its walls still carry fringe to it, which pulls it down."""
    from test_equilibrium import plate_on_beams

    plate = plate_on_beams()
    model = plate.capacitance(ElectrostaticParams(interaction_distance=30 * UM))
    over = plate.electrode_at((50 * UM, 50 * UM), "device")
    under = plate.electrode_at((50 * UM, 50 * UM), "handle")
    assert model.clusters == ()
    assert -model.fringe[over.index, under.index] > 0
    assert -model.fringe[over.index, under.index] < -model.released[over.index, under.index]
    up = model.dofs.index((over.body, "uz"))
    assert model.fringe[over.index, under.index] != 0
    assert model.gradient[over.index, under.index, up] != 0
    alone = plate.capacitance(
        ElectrostaticParams(interaction_distance=30 * UM, include_fringe=False)
    )
    assert not alone.fringe.any()


def facing_pads() -> Structure:
    """Two released pads 40 um apart, each on a beam to its own anchor, hanging over a
    gap in the handle so that neither has handle beneath; they are in no cluster at
    an interaction distance of 10 um."""
    device = [
        box(-500, -100, -300, 100),
        box(-300, -5, -100, 5),
        box(-100, -50, 0, 50),
        box(40, -50, 140, 50),
        box(140, -5, 340, 5),
        box(340, -100, 540, 100),
    ]
    handle = Polygon(box(-600, -200, 600, 200).exterior, [box(-280, -150, 320, 150).exterior])
    layout = Layout.from_shapely(device=device, handle=[handle], unit=UM)
    stack = Stack(
        device=Layer(thickness=THICKNESS, material=Material.silicon()),
        box=Layer(thickness=BOX, material=Material.oxide()),
        handle=Layer(thickness=400 * UM, material=Material.silicon()),
    )
    undercut = ExtractionParams(
        release_hole_max_area=100 * UM * UM, release_distance=7 * UM, body_min_width=20 * UM
    )
    return extract(layout, stack, undercut)


def test_bodies_of_different_clusters_closing_on_each_other_are_seen() -> None:
    """Electrodes farther apart than the interaction distance at rest are in no cluster,
    yet their gap counts in the clearance once they come within reach, and touching
    raises as within a cluster."""
    pads = facing_pads()
    model = pads.capacitance(ElectrostaticParams(interaction_distance=10 * UM))
    left = pads.electrode_at((-50 * UM, 0.0), "device")
    right = pads.electrode_at((90 * UM, 0.0), "device")
    assert model.clusters == ()
    assert model.clearance == 1.0
    toward = model.dofs.index((left.body, "ux"))
    place = np.zeros(len(model.dofs))
    place[toward] = 15 * UM
    far = model.at(place)
    assert far.clearance == 1.0  # still beyond reach
    place[toward] = 32 * UM
    near = model.at(place)
    assert near.clearance == pytest.approx(8 / 40)
    assert near.closest == tuple(sorted((left.index, right.index)))
    place[toward] = 41 * UM
    with pytest.raises(ValueError, match="touch"):
        model.at(place)


def test_turn_that_swings_finger_tips_into_the_stator_is_seen(
    comb: Structure, at: dict[str, int]
) -> None:
    """Clearance and touching follow a turn of the shuttle as they follow its travel: the
    tips of its fingers, far from its center, sweep across the gap."""
    rest = comb.capacitance(
        ElectrostaticParams(
            dense_cutoff=0.0, plates_cutoff=0.0, interaction_distance=10 * UM, include_fringe=False
        )
    )
    turned = rest.at([0.0, 0.0, 0.0, 0.0, 0.0, 0.004])
    assert 0.0 < turned.clearance < 1.0
    assert turned.closest == (at["stator"], at["shuttle"])
    with pytest.raises(ValueError, match="touch"):
        rest.at([0.0, 0.0, 0.0, 0.0, 0.0, 0.05])


def two_combs() -> Structure:
    """Two copies of the comb one above the other, each shuttle its own body."""
    device = []
    for shift in (0.0, 600.0):
        stator = [
            box(300, shift + 100 + PITCH * k, 360, shift + 100 + PITCH * k + WIDTH)
            for k in range(FINGERS + 1)
        ]
        rotor = [
            box(
                320,
                shift + 100 + PITCH * k + WIDTH + GAP,
                380,
                shift + 100 + PITCH * k + 2 * WIDTH + GAP,
            )
            for k in range(FINGERS)
        ]
        spine = box(380, shift + 60, 420, shift + 200)
        beams = [box(420, shift + 70, 620, shift + 80), box(420, shift + 180, 620, shift + 190)]
        device += [box(0, shift, 300, shift + 300), *stator, *rotor, spine, *beams]
        device.append(box(620, shift, 920, shift + 300))
    layout = Layout.from_shapely(device=device, handle=[box(-50, -50, 970, 950)], unit=UM)
    stack = Stack(
        device=Layer(thickness=THICKNESS, material=Material.silicon()),
        box=Layer(thickness=BOX, material=Material.oxide()),
        handle=Layer(thickness=400 * UM, material=Material.silicon()),
    )
    params = ExtractionParams(
        release_hole_max_area=100 * UM * UM, release_distance=60 * UM, body_min_width=20 * UM
    )
    return extract(layout, stack, params)


def test_clusters_whose_bodies_did_not_move_keep_their_rest_solution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An exact evaluation solves the field of a cluster only where one of its bodies
    has moved in the plane; the others keep the solution found at rest."""
    from connflex_core.electrostatics import capacitance as module
    from connflex_core.electrostatics.sidewalls import solve_sidewalls as original

    structure = two_combs()
    params = ElectrostaticParams(
        element_size=4 * UM, interaction_distance=10 * UM, include_fringe=False
    )
    model = structure.capacitance(params)
    assert len(model.clusters) == 2
    solved: list[int] = []

    def counted(*args: Any, **kwargs: Any) -> Any:
        solved.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(module, "solve_sidewalls", counted)
    lower = next(k for k, (body, name) in enumerate(model.dofs) if name == "ux")
    place = np.zeros(len(model.dofs))
    place[lower] = 0.5 * UM
    moved = model.at(place, exact=True)
    assert solved == [1]
    assert moved.at(np.zeros(len(model.dofs)), exact=True).sidewall == pytest.approx(
        model.sidewall, rel=1e-12, abs=1e-30
    )
    assert solved == [1]
