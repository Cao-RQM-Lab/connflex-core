"""Threads that run independent pieces of an analysis."""

from __future__ import annotations

import os
import threading
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Generic, TypeVar

from connflex_core._validation import require_at_least

T = TypeVar("T")
K = TypeVar("K")

_here = threading.local()


class Workers:
    """``count`` threads for pieces of work that do not depend on one another.

    None takes every processor of the machine. With one worker, with one piece, or from
    inside a worker, pieces run one after another in the calling thread.
    """

    def __init__(self, count: int | None = None) -> None:
        self.count = (os.cpu_count() or 1) if count is None else count
        require_at_least("workers", self.count, 1)

    def run(
        self, pieces: Sequence[Callable[[], T]], done: Callable[[], None] | None = None
    ) -> list[T]:
        """Results of ``pieces`` in their order; ``done`` is called in the calling thread
        as each piece finishes."""
        if self.count == 1 or len(pieces) < 2 or getattr(_here, "inside", False):
            results = []
            for piece in pieces:
                results.append(piece())
                if done is not None:
                    done()
            return results
        with ThreadPoolExecutor(min(self.count, len(pieces))) as pool:
            futures = [pool.submit(_inside, piece) for piece in pieces]
            for _ in as_completed(futures):
                if done is not None:
                    done()
            return [future.result() for future in futures]


def _inside(piece: Callable[[], T]) -> T:
    _here.inside = True
    return piece()


class Memo(Generic[K, T]):
    """Values computed once per key: the first thread to ask computes, the others wait."""

    def __init__(self) -> None:
        self.values: dict[K, T] = {}
        self._pending: dict[K, threading.Event] = {}
        self._lock = threading.Lock()

    def get(self, key: K, compute: Callable[[], T]) -> T:
        """Value of ``key``, from ``compute`` when no thread has found it yet."""
        while True:
            with self._lock:
                if key in self.values:
                    return self.values[key]
                waiting = self._pending.get(key)
                mine = waiting is None
                if waiting is None:
                    waiting = self._pending[key] = threading.Event()
            if not mine:
                waiting.wait()
                continue
            try:
                value = compute()
            except BaseException:
                with self._lock:
                    del self._pending[key]
                waiting.set()
                raise
            with self._lock:
                self.values[key] = value
                del self._pending[key]
            waiting.set()
            return value
