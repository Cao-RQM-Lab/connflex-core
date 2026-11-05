import json
from dataclasses import replace

import numpy as np
import pytest
from shapely import Polygon, affinity, box

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
    return ElectrostaticParams(
        sidewall_method="field", element_size=4 * UM, interaction_distance=30 * UM
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

    with pytest.raises(ValueError, match="sidewall_method"):
        comb.capacitance(ElectrostaticParams(element_size=4 * UM, interaction_distance=30 * UM))
    with pytest.raises(ValueError, match="interaction_distance"):
        comb.capacitance(replace(settings, interaction_distance=None))
    with pytest.raises(ValueError, match="element_size"):
        comb.capacitance(replace(settings, element_size=None))

    assert comb.capacitance(replace(settings, interaction_distance=1 * UM)).clusters == ()


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
        sidewall_method="plates", interaction_distance=10 * UM, include_plates=False
    )
    model = comb.capacitance(params)
    stator, shuttle = at["stator"], at["shuttle"]
    assert -model.sidewall[stator, shuttle] == pytest.approx(
        GAPS * EPS0 * THICKNESS * OVERLAP / GAP, rel=1e-9
    )
    (cluster,) = model.clusters
    assert (cluster.electrodes, cluster.elements, cluster.windows) == ((stator, shuttle), 16, 0)
    assert cluster.walls.gap == pytest.approx(GAP * UM)

    potentials = np.zeros(len(comb.electrodes))
    potentials[stator] = 10.0
    force, stiffness = model.force(potentials), model.stiffness(potentials)
    drive = GAPS * 0.5 * EPS0 * THICKNESS * 10.0**2 / (GAP * UM)
    assert force[:2] == pytest.approx([-drive, 0.0], rel=1e-9, abs=1e-9 * drive)
    side = GAPS * EPS0 * THICKNESS * OVERLAP * UM * 10.0**2 / (GAP * UM) ** 3
    assert stiffness[1, 1] == pytest.approx(-side, rel=1e-9)


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
            sidewall_method="field", element_size=4 * UM, interaction_distance=30 * UM
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
        }
    ]
    assert report["unassigned_area"] == pytest.approx(2 * 200 * 10 * UM**2)
