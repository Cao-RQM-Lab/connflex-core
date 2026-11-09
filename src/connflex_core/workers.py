"""Threads and processes that run independent pieces of an analysis."""

from __future__ import annotations

import atexit
import contextlib
import os
import pickle
import signal
import struct
import subprocess
import sys
import threading
import traceback
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from typing import IO, Any, Generic, TypeVar

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
        pool = ThreadPoolExecutor(min(self.count, len(pieces)))
        try:
            futures = [pool.submit(_inside, piece) for piece in pieces]
            for _ in as_completed(futures):
                if done is not None:
                    done()
            return [future.result() for future in futures]
        finally:  # after an interrupt, pieces not yet started stay so
            pool.shutdown(wait=False, cancel_futures=True)

    def apply(self, function: Callable[..., T], arguments: tuple[Any, ...]) -> T:
        """``function(*arguments)`` in a child process, or in this one with one worker.

        ``function`` must be importable by name, so defined in a module and not in
        ``__main__``, and its arguments and result picklable. An exception of the
        function is raised here with its type, a note holding its traceback in the
        child; a child that dies raises ``RuntimeError``.
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


HEADER = struct.Struct("<Q")
LONGEST = 1 << 40  # a message longer than this announces a stream out of step


def _pack(message: object) -> bytes:
    return pickle.dumps(message, protocol=pickle.HIGHEST_PROTOCOL)


def _write(stream: IO[bytes], data: bytes) -> None:
    """Send one message: its length, then its bytes, in one write."""
    stream.write(HEADER.pack(len(data)) + data)
    stream.flush()


def _read(stream: IO[bytes]) -> bytes:
    """Bytes of the next message; ``EOFError`` when the stream has ended."""
    header = stream.read(HEADER.size)
    if len(header) < HEADER.size:
        raise EOFError
    (length,) = HEADER.unpack(header)
    if length > LONGEST:
        raise pickle.UnpicklingError(f"a message of {length} bytes is announced")
    data = stream.read(length)
    if len(data) < length:
        raise EOFError
    return data


def _receive(stream: IO[bytes]) -> Any:
    return pickle.loads(_read(stream))


class _Child:
    """One child interpreter running ``serve``, spoken to over pipes in messages of a
    length then a pickle, so that a message either side cannot read costs it nothing
    but that message."""

    def __init__(self) -> None:
        self.process = subprocess.Popen(
            [sys.executable, "-m", "connflex_core._serve"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        assert self.process.stdin is not None
        assert self.process.stdout is not None
        self.outbox: IO[bytes] = self.process.stdin
        self.inbox: IO[bytes] = self.process.stdout
        _write(self.outbox, _pack(sys.path))

    def ask(self, message: bytes) -> tuple[bool, Any]:
        """Reply of the child to a packed ``(function, arguments)``: whether the function
        returned, and its result or its exception. ``EOFError`` when the child has died."""
        try:
            _write(self.outbox, message)
            reply: tuple[bool, Any] = _receive(self.inbox)
        except (OSError, EOFError, pickle.UnpicklingError) as exc:
            raise EOFError from exc
        return reply

    def alive(self) -> bool:
        return self.process.poll() is None

    def close(self) -> None:
        """Let the child finish what it has and leave."""
        try:
            self.outbox.close()
            self.process.wait(timeout=5.0)
        except (OSError, subprocess.TimeoutExpired):
            self.process.kill()

    def kill(self) -> None:
        self.process.kill()
        self.process.wait()


class _Pool:
    """Idle and busy children of this interpreter, started as they are first needed."""

    def __init__(self) -> None:
        self.idle: list[_Child] = []
        self.alive = 0
        self.room = threading.Condition()

    def apply(self, function: Callable[..., T], arguments: tuple[Any, ...], limit: int) -> T:
        """``function(*arguments)`` in an idle child, or a new one while fewer than
        ``limit`` are alive; otherwise wait for one."""
        if getattr(function, "__module__", None) == "__main__":
            raise ValueError(
                f"{function.__name__} is defined in __main__, which a worker process cannot "
                f"import: define it in a module"
            )
        message = _pack((function, arguments))
        child = self._take(limit)
        try:
            succeeded, value = child.ask(message)
        except EOFError:
            self._drop(child)
            raise RuntimeError(
                "a worker process died or lost its place; its error, if any, is above"
            ) from None
        except BaseException:
            self._drop(child)  # its reply, if one comes, would answer the next question
            raise
        with self.room:
            self.idle.append(child)
            self.room.notify()
        if succeeded:
            result: T = value
            return result
        raise value

    def _take(self, limit: int) -> _Child:
        """An idle child that is alive, or a new one while fewer than ``limit`` are alive."""
        with self.room:
            while True:
                while self.idle:
                    child = self.idle.pop()
                    if child.alive():
                        return child
                    self.alive -= 1  # died while idle
                if self.alive < limit:
                    self.alive += 1
                    break
                self.room.wait()
        try:
            return _Child()
        except BaseException:
            with self.room:
                self.alive -= 1
                self.room.notify()
            raise

    def _drop(self, child: _Child) -> None:
        child.kill()
        with self.room:
            self.alive -= 1
            self.room.notify()

    def close(self) -> None:
        with self.room:
            children, self.idle = self.idle, []
            self.alive -= len(children)
            self.room.notify_all()
        for child in children:
            child.close()

    def forget(self) -> None:
        """Start afresh in a process forked from this one: its children are the parent's."""
        self.idle, self.alive, self.room = [], 0, threading.Condition()


_POOL = _Pool()
atexit.register(_POOL.close)
if hasattr(os, "register_at_fork"):
    os.register_at_fork(after_in_child=_POOL.forget)


def serve() -> None:
    """Run pieces that a parent sends over stdin until it closes the pipe: the child's
    half of ``Workers``, which ``python -m connflex_core._serve`` starts."""
    inbox: IO[bytes] = sys.stdin.buffer
    outbox: IO[bytes] = os.fdopen(os.dup(1), "wb")
    os.dup2(2, 1)  # whatever a piece writes to standard output, even from C, goes to standard error
    sys.stdout = sys.stderr
    signal.signal(signal.SIGINT, signal.SIG_IGN)  # the parent handles interrupts of its terminal
    sys.path[:] = _receive(inbox)
    while True:
        try:
            data = _read(inbox)
        except EOFError:
            return
        try:
            function, arguments = pickle.loads(data)
        except Exception as exc:
            _write(outbox, _pack((False, RuntimeError(f"a worker could not read a piece: {exc}"))))
            continue
        try:
            message: tuple[bool, object] = (True, function(*arguments))
        except BaseException as exc:
            with contextlib.suppress(Exception):
                exc.add_note(
                    f"Traceback in worker process {os.getpid()} (most recent call last):\n"
                    + "".join(traceback.format_tb(exc.__traceback__)).rstrip()
                )
            message = (False, exc)
        try:
            packed = _pack(message)
        except Exception as exc:
            failed = RuntimeError(f"result of {function.__name__} could not be sent: {exc}")
            packed = _pack((False, failed))
        _write(outbox, packed)
