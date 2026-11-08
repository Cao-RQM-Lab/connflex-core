"""Threads and processes that run independent pieces of an analysis."""

from __future__ import annotations

import atexit
import os
import pickle
import subprocess
import sys
import threading
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from typing import Any, Generic, TypeVar

from connflex_core._validation import require_at_least

T = TypeVar("T")
K = TypeVar("K")

_here = threading.local()


class Workers:
    """``count`` threads and processes for pieces of work that do not depend on one another.

    None takes every processor of the machine. ``run`` runs pieces in threads, for work
    that releases the interpreter's lock or is short; ``solve`` and ``apply`` run a
    module-level function in child interpreters running ``serve``, for work that holds
    it, sharing one pool of children that the first call starts and the interpreter
    keeps until it exits. With one worker everything runs in the calling thread and process.
    """

    def __init__(self, count: int | None = None) -> None:
        self.count = (os.cpu_count() or 1) if count is None else count
        require_at_least("workers", self.count, 1)

    def run(
        self, pieces: Sequence[Callable[[], T]], done: Callable[[], None] | None = None
    ) -> list[T]:
        """Results of ``pieces`` in their order; ``done`` is called in the calling thread
        as each piece finishes. With one worker, with one piece, or from inside a worker,
        pieces run one after another in the calling thread."""
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

    def apply(self, function: Callable[..., T], arguments: tuple[Any, ...]) -> T:
        """``function(*arguments)`` in a child process, or in this one with one worker.

        ``function`` must be importable by name and its arguments and result picklable.
        An exception of the function is raised here with its type; a child that dies
        raises ``RuntimeError``.
        """
        if self.count == 1:
            return function(*arguments)
        return _POOL.apply(function, arguments, self.count)

    def solve(
        self,
        function: Callable[..., T],
        arguments: Sequence[tuple[Any, ...]],
        done: Callable[[], None] | None = None,
    ) -> list[T]:
        """Results of ``function`` on each tuple of ``arguments`` in their order, through
        ``apply`` by ``run``; in this process when there is one worker or one tuple."""
        if self.count == 1 or len(arguments) < 2:
            results = []
            for tuple_ in arguments:
                results.append(function(*tuple_))
                if done is not None:
                    done()
            return results
        return self.run([partial(self.apply, function, tuple_) for tuple_ in arguments], done)


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


class _Child:
    """One child interpreter running ``serve``, spoken to over pipes."""

    def __init__(self) -> None:
        self.process = subprocess.Popen(
            [sys.executable, "-m", "connflex_core._serve"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        self._send(sys.path)

    def _send(self, message: object) -> None:
        assert self.process.stdin is not None
        pickle.dump(message, self.process.stdin, protocol=pickle.HIGHEST_PROTOCOL)
        self.process.stdin.flush()

    def ask(self, function: Callable[..., T], arguments: tuple[Any, ...]) -> T:
        """Run ``function(*arguments)`` in the child; ``EOFError`` when it has died."""
        assert self.process.stdout is not None
        try:
            self._send((function, arguments))
            succeeded, value = pickle.load(self.process.stdout)
        except (OSError, EOFError, pickle.UnpicklingError) as exc:
            raise EOFError from exc
        if succeeded:
            result: T = value
            return result
        raise value

    def close(self) -> None:
        assert self.process.stdin is not None
        try:
            self.process.stdin.close()
            self.process.wait(timeout=5.0)
        except (OSError, subprocess.TimeoutExpired):
            self.process.kill()


class _Pool:
    """Idle and busy children of this interpreter, started as they are first needed."""

    def __init__(self) -> None:
        self.idle: list[_Child] = []
        self.alive = 0
        self.room = threading.Condition()

    def apply(self, function: Callable[..., T], arguments: tuple[Any, ...], limit: int) -> T:
        """``function(*arguments)`` in an idle child, or a new one while fewer than
        ``limit`` are alive; otherwise wait for one."""
        with self.room:
            while not self.idle and self.alive >= limit:
                self.room.wait()
            child = self.idle.pop() if self.idle else None
            if child is None:
                self.alive += 1
        if child is None:
            try:
                child = _Child()
            except BaseException:
                with self.room:
                    self.alive -= 1
                    self.room.notify()
                raise
        healthy = True
        try:
            return child.ask(function, arguments)
        except EOFError:
            healthy = False
            child.close()
            raise RuntimeError("a worker process died; its error, if any, is above") from None
        finally:
            with self.room:
                if healthy:
                    self.idle.append(child)
                else:
                    self.alive -= 1
                self.room.notify()

    def close(self) -> None:
        with self.room:
            children, self.idle = self.idle, []
        for child in children:
            child.close()


_POOL = _Pool()
atexit.register(_POOL.close)


def serve() -> None:
    """Run pieces that a parent sends over stdin until it closes the pipe: the child's
    half of ``Workers``, which ``python -m connflex_core._serve`` starts."""
    inbox, outbox = sys.stdin.buffer, sys.stdout.buffer
    sys.stdout = sys.stderr  # nothing a piece prints may reach the pipe
    sys.path[:] = pickle.load(inbox)
    while True:
        try:
            function, arguments = pickle.load(inbox)
        except EOFError:
            return
        try:
            message: tuple[bool, object] = (True, function(*arguments))
        except BaseException as exc:
            message = (False, exc)
        try:
            pickle.dump(message, outbox, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception as exc:
            failed = RuntimeError(f"result of {function.__name__} could not be sent: {exc}")
            pickle.dump((False, failed), outbox, protocol=pickle.HIGHEST_PROTOCOL)
        outbox.flush()
