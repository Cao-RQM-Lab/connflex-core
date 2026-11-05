import logging
import re
from collections.abc import Callable
from dataclasses import dataclass, replace

import pytest
from shapely import box

from connflex_core import (
    ElectrostaticParams,
    ExtractionParams,
    Layout,
    Stack,
    extract,
    sidewall_capacitance,
)
from connflex_core.progress import stage

EXTRACT_STAGES = [
    "heal",
    "release",
    "prune",
    "nets",
    "opening",
    "bodies",
    "candidates",
    "beams",
    "mass",
    "electrodes",
]
PAD_A, BEAM, PAD_B = box(0, 0, 100, 100), box(100, 45, 300, 55), box(300, 0, 400, 100)
HANDLE_A = box(-50, -50, 150, 150)


@dataclass
class Run:
    name: str
    total: int | None
    done: int = 0
    finished: bool = False


class Recorder:
    def __init__(self) -> None:
        self.runs: list[Run] = []

    def start(self, stage: str, total: int | None) -> None:
        assert not self.runs or self.runs[-1].finished
        self.runs.append(Run(stage, total))

    def advance(self, n: int = 1) -> None:
        assert not self.runs[-1].finished
        self.runs[-1].done += n

    def finish(self) -> None:
        assert not self.runs[-1].finished
        self.runs[-1].finished = True


def test_package_logger_is_silent_without_a_handler() -> None:
    logger = logging.getLogger("connflex_core")
    assert any(isinstance(handler, logging.NullHandler) for handler in logger.handlers)


def test_stage_logs_both_ends_and_forwards_events(caplog: pytest.LogCaptureFixture) -> None:
    recorder = Recorder()
    with (
        caplog.at_level(logging.INFO, logger="connflex_core"),
        stage("heal", recorder, total=3) as current,
    ):
        current.advance(2)
        current.advance()
        current.count(device=5, handle=1)
    assert recorder.runs == [Run("heal", 3, 3, True)]
    messages = [record.getMessage() for record in caplog.records]
    assert messages[0] == "heal: started"
    assert re.fullmatch(r"heal: done in \d+\.\d+ s, device 5, handle 1", messages[1])


def test_stage_finishes_the_callback_when_the_body_raises(
    caplog: pytest.LogCaptureFixture,
) -> None:
    recorder = Recorder()
    with (
        caplog.at_level(logging.INFO, logger="connflex_core"),
        pytest.raises(RuntimeError),
        stage("release", recorder),
    ):
        raise RuntimeError
    assert recorder.runs == [Run("release", None, 0, True)]
    assert [record.getMessage() for record in caplog.records] == ["release: started"]


def test_stage_without_a_callback_only_logs(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.INFO, logger="connflex_core"), stage("prune", None) as current:
        current.advance()
    assert len(caplog.records) == 2


def test_extract_reports_every_stage_in_order(
    stack: Stack,
    params: Callable[..., ExtractionParams],
    layout_um: Callable[..., Layout],
    caplog: pytest.LogCaptureFixture,
) -> None:
    recorder = Recorder()
    with caplog.at_level(logging.INFO, logger="connflex_core"):
        structure = extract(
            layout_um([PAD_A, BEAM, PAD_B], [HANDLE_A]), stack, params(), progress=recorder
        )
    assert [run.name for run in recorder.runs] == EXTRACT_STAGES
    assert all(run.finished for run in recorder.runs)
    assert all(run.total is None or run.done == run.total for run in recorder.runs)
    totals = {run.name: run.total for run in recorder.runs}
    assert totals["bodies"] == sum(
        len(body.device.geoms) + len(body.box.geoms) + len(body.handle.geoms)
        for body in structure.bodies
    )
    assert totals["mass"] == len(structure.bodies)
    assert totals["beams"] == len(structure.beams) + len(structure.rejected)
    assert len(caplog.records) == 2 * len(EXTRACT_STAGES)


def test_modal_and_static_report_assemble_and_solve(
    stack: Stack, params: Callable[..., ExtractionParams], layout_um: Callable[..., Layout]
) -> None:
    structure = extract(layout_um([PAD_A, BEAM, PAD_B], [HANDLE_A]), stack, params())
    for solve in (
        lambda recorder: structure.modal(progress=recorder),
        lambda recorder: structure.static({1: [1e-6, 0, 0, 0, 0, 0]}, progress=recorder),
    ):
        recorder = Recorder()
        solve(recorder)
        assert recorder.runs == [Run("assemble", None, 0, True), Run("solve", None, 0, True)]


def test_plates_report_one_stage_that_advances_once_per_device_electrode(
    stack: Stack, params: Callable[..., ExtractionParams], layout_um: Callable[..., Layout]
) -> None:
    structure = extract(layout_um([PAD_A, BEAM, PAD_B], [HANDLE_A]), stack, params())
    recorder = Recorder()
    assert len(structure.plates(progress=recorder)) == 1
    assert recorder.runs == [Run("plates", 2, 2, True)]


def test_sidewall_capacitance_reports_its_stages() -> None:
    recorder = Recorder()
    squares = [box(0, 0, 10e-6, 10e-6), box(12e-6, 0, 22e-6, 10e-6)]
    field = ElectrostaticParams(
        element_size=1e-6,
        interaction_distance=5e-6,
        dense_cutoff=float("inf"),
        plates_cutoff=float("inf"),
    )
    result = sidewall_capacitance(squares, field, progress=recorder)
    elements = len(result.mesh.conductor)
    assert recorder.runs == [
        Run("facing", None, 0, True),
        Run("mesh", None, 0, True),
        Run("assemble", elements, elements, True),
        Run("solve", None, 0, True),
        Run("derivatives", 1, 1, True),
        Run("curvature", None, 0, True),
    ]
    recorder = Recorder()
    limited = replace(field, max_elements=elements - 1, window_elements=elements // 4)
    windows = sidewall_capacitance(squares, limited, progress=recorder).windows
    assert recorder.runs == [
        Run("facing", None, 0, True),
        Run("mesh", None, 0, True),
        Run("windows", windows, windows, True),
    ]
    recorder = Recorder()
    plates = ElectrostaticParams(dense_cutoff=0.0, plates_cutoff=0.0, interaction_distance=5e-6)
    sidewall_capacitance(squares, plates, progress=recorder)
    assert recorder.runs == [Run("facing", None, 0, True)]
