import logging
import re

import pytest

from connflex_core.progress import stage


class Recorder:
    def __init__(self) -> None:
        self.events: list[tuple[object, ...]] = []

    def start(self, stage: str, total: int | None) -> None:
        self.events.append(("start", stage, total))

    def advance(self, n: int = 1) -> None:
        self.events.append(("advance", n))

    def finish(self) -> None:
        self.events.append(("finish",))


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
    assert recorder.events == [("start", "heal", 3), ("advance", 2), ("advance", 1), ("finish",)]
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
    assert recorder.events == [("start", "release", None), ("finish",)]
    assert [record.getMessage() for record in caplog.records] == ["release: started"]


def test_stage_without_a_callback_only_logs(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.INFO, logger="connflex_core"), stage("prune", None) as current:
        current.advance()
    assert len(caplog.records) == 2
