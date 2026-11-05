"""Stage logging and progress reporting."""

from __future__ import annotations

import logging
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Protocol

logger = logging.getLogger("connflex_core")
logger.addHandler(logging.NullHandler())


class Progress(Protocol):
    """Receiver of stage events from ``extract`` and from every analysis.

    ``start`` opens a stage with the number of ``advance`` steps to expect, or None when
    the stage is one indivisible operation. ``finish`` closes it, also when the stage fails.
    """

    def start(self, stage: str, total: int | None) -> None: ...

    def advance(self, n: int = 1) -> None: ...

    def finish(self) -> None: ...


class _Silent:
    def start(self, stage: str, total: int | None) -> None:
        pass

    def advance(self, n: int = 1) -> None:
        pass

    def finish(self) -> None:
        pass


@dataclass
class Stage:
    """One running stage: forwards ``advance`` and collects counts for the finish log line."""

    name: str
    progress: Progress
    counts: dict[str, int] = field(default_factory=dict)

    def advance(self, n: int = 1) -> None:
        self.progress.advance(n)

    def count(self, **counts: int) -> None:
        self.counts.update(counts)


@contextmanager
def stage(name: str, progress: Progress | None, total: int | None = None) -> Iterator[Stage]:
    """Run a stage: log start and finish with elapsed time, and report to ``progress``."""
    current = Stage(name, _Silent() if progress is None else progress)
    logger.info("%s: started", name)
    current.progress.start(name, total)
    started = time.perf_counter()
    try:
        yield current
    finally:
        current.progress.finish()
    elapsed = time.perf_counter() - started
    details = "".join(f", {key} {value}" for key, value in current.counts.items())
    logger.info("%s: done in %.2f s%s", name, elapsed, details)


@contextmanager
def quiet(name: str, total: int | None = None) -> Iterator[Stage]:
    """Run a stage inside another one: nothing is logged and nothing is reported."""
    yield Stage(name, _Silent())
