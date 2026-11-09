import os
import pickle
import signal
import sys
import threading
import time
import types
from collections.abc import Callable
from functools import partial
from typing import Any

import pytest

from connflex_core.workers import Memo, Workers


def _slow(k: int) -> int:
    time.sleep(0.01 * (3 - k))
    return k


def test_results_keep_the_order_of_pieces_and_each_piece_is_reported() -> None:
    reported: list[int] = []
    pieces = [partial(_slow, k) for k in range(3)]
    assert Workers(3).run(pieces, lambda: reported.append(1)) == [0, 1, 2]
    assert len(reported) == 3
    assert Workers(1).run(pieces, lambda: reported.append(1)) == [0, 1, 2]
    assert len(reported) == 6


def test_work_from_inside_a_worker_runs_in_that_thread() -> None:
    workers = Workers(2)

    def nested() -> bool:
        mine = threading.get_ident()
        return all(ident == mine for ident in workers.run([threading.get_ident] * 4))

    assert workers.run([nested, nested]) == [True, True]


def test_error_of_a_piece_reaches_the_caller() -> None:
    def fail() -> int:
        raise ValueError("piece")

    with pytest.raises(ValueError, match="piece"):
        Workers(2).run([lambda: 1, fail])


def test_workers_must_be_at_least_one() -> None:
    with pytest.raises(ValueError, match="workers"):
        Workers(0)
    assert Workers(None).count >= 1


def test_memo_computes_each_key_once_for_concurrent_callers() -> None:
    memo: Memo[str, int] = Memo()
    computed: list[int] = []
    barrier = threading.Barrier(4)

    def ask() -> int:
        barrier.wait()

        def compute() -> int:
            computed.append(1)
            time.sleep(0.02)
            return 7

        return memo.get("key", compute)

    assert Workers(4).run([ask] * 4) == [7] * 4
    assert computed == [1]
    assert memo.get("key", lambda: 0) == 7


def _where(k: int) -> tuple[int, int]:
    return os.getpid(), k * k


def _fail_on_one(k: int) -> int:
    if k == 1:
        raise ValueError(f"piece {k}")
    return k


def _quit(k: int) -> int:
    os._exit(3)


def test_solve_runs_a_module_function_in_other_processes_in_order() -> None:
    reported: list[int] = []
    results = Workers(2).solve(_where, [(k,) for k in range(4)], lambda: reported.append(1))
    assert [square for _, square in results] == [0, 1, 4, 9]
    assert all(pid != os.getpid() for pid, _ in results)
    assert len(reported) == 4


def test_solve_stays_in_this_process_with_one_worker_or_one_piece() -> None:
    assert Workers(1).solve(_where, [(2,), (3,)]) == [(os.getpid(), 4), (os.getpid(), 9)]
    assert Workers(3).solve(_where, [(2,)]) == [(os.getpid(), 4)]


def _where_after_a_moment(k: int) -> tuple[int, int]:
    time.sleep(0.05)  # long enough for both pieces to hold a process at once
    return _where(k)


def test_processes_are_kept_for_later_work() -> None:
    first = {pid for pid, _ in Workers(2).solve(_where_after_a_moment, [(0,), (1,)])}
    second = {pid for pid, _ in Workers(2).solve(_where_after_a_moment, [(0,), (1,)])}
    assert len(first) == 2
    assert second == first


def test_error_in_a_process_reaches_the_caller_with_its_type() -> None:
    with pytest.raises(ValueError, match="piece 1"):
        Workers(2).solve(_fail_on_one, [(0,), (1,)])
    assert Workers(2).solve(_fail_on_one, [(0,), (2,)]) == [0, 2]


def test_a_process_that_dies_is_reported_and_replaced() -> None:
    with pytest.raises(RuntimeError, match="worker process"):
        Workers(2).solve(_quit, [(0,), (1,)])
    assert [square for _, square in Workers(2).solve(_where, [(2,), (3,)])] == [4, 9]


def test_solve_from_inside_a_thread_worker_still_uses_processes() -> None:
    workers = Workers(2)

    def nested() -> list[tuple[int, int]]:
        return workers.solve(_where, [(1,), (2,)])

    for found in workers.run([nested, nested]):
        assert [square for _, square in found] == [1, 4]
        assert all(pid != os.getpid() for pid, _ in found)


class _Interrupted(BaseException):
    """Stands in for an interrupt of the calling thread while it waits on a child."""


def test_interrupt_while_waiting_discards_the_child_and_later_results_are_right(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A child whose reply was never read is out of step with this process: it is
    killed, not kept, so that no later piece reads the reply of an earlier one."""
    from connflex_core import workers as module

    workers = Workers(2)
    workers.apply(_where, (0,))
    receive = module._receive

    def interrupted(stream: Any) -> Any:
        monkeypatch.setattr(module, "_receive", receive)
        raise _Interrupted

    monkeypatch.setattr(module, "_receive", interrupted)
    before = module._POOL.alive
    with pytest.raises(_Interrupted):
        workers.apply(_where_after_a_moment, (1,))
    assert module._POOL.alive == before - 1
    assert workers.apply(_where, (7,))[1] == 49
    assert workers.apply(_where, (8,))[1] == 64


def test_a_piece_that_cannot_be_sent_leaves_the_child_in_step() -> None:
    from connflex_core.workers import _POOL

    workers = Workers(2)
    workers.apply(_where, (0,))
    idle = len(_POOL.idle)
    with pytest.raises((pickle.PicklingError, AttributeError)):
        workers.apply(_where, (lambda: 0,))
    assert len(_POOL.idle) == idle
    assert workers.apply(_where, (3,))[1] == 9


def _unsendable(k: int) -> Callable[[], int]:
    return lambda: k


def test_a_result_that_cannot_be_sent_is_an_error_and_the_child_goes_on() -> None:
    workers = Workers(2)
    with pytest.raises(RuntimeError, match="could not be sent"):
        workers.apply(_unsendable, (1,))
    assert workers.apply(_where, (3,))[1] == 9


def test_a_function_defined_in_main_is_refused_before_it_is_sent() -> None:
    def mine(k: int) -> int:
        return k

    mine.__module__ = "__main__"
    with pytest.raises(ValueError, match="__main__"):
        Workers(2).apply(mine, (1,))


def test_a_piece_the_child_cannot_read_is_an_error_and_the_child_goes_on() -> None:
    """A class that only this process knows pickles by reference and cannot be read in
    the child, which says so and stays in step for the next piece."""
    module = types.ModuleType("only_here")
    module.Thing = type("Thing", (), {"__module__": "only_here"})  # type: ignore[attr-defined]
    sys.modules["only_here"] = module
    try:
        workers = Workers(2)
        pid, _ = workers.apply(_where, (0,))
        with pytest.raises(RuntimeError, match="could not read"):
            workers.apply(_where, (module.Thing(),))
        assert workers.apply(_where, (3,)) == (pid, 9)
    finally:
        del sys.modules["only_here"]


def test_close_forgets_its_children_and_the_pool_starts_anew() -> None:
    from connflex_core.workers import _POOL

    workers = Workers(2)
    workers.solve(_where_after_a_moment, [(0,), (1,)])
    _POOL.close()
    assert (_POOL.alive, _POOL.idle) == (0, [])
    assert workers.apply(_where, (2,))[1] == 4


def test_children_killed_while_idle_are_replaced() -> None:
    from connflex_core.workers import _POOL

    workers = Workers(2)
    workers.solve(_where_after_a_moment, [(0,), (1,)])
    for child in _POOL.idle:
        child.kill()
    assert [square for _, square in workers.solve(_where, [(2,), (3,)])] == [4, 9]


@pytest.mark.skipif(sys.platform == "win32", reason="signals of a terminal are POSIX")
def test_an_interrupt_from_the_terminal_leaves_idle_children_alive() -> None:
    """A terminal sends its interrupt to the parent and its children alike; the
    children leave it to the parent and stay for the next call."""
    from connflex_core.workers import _POOL

    workers = Workers(2)
    first = {pid for pid, _ in workers.solve(_where_after_a_moment, [(0,), (1,)])}
    for child in _POOL.idle:
        os.kill(child.process.pid, signal.SIGINT)
    time.sleep(0.2)
    assert all(child.alive() for child in _POOL.idle)
    second = {pid for pid, _ in workers.solve(_where_after_a_moment, [(0,), (1,)])}
    assert second == first


def test_interrupt_of_the_calling_thread_cancels_pieces_not_yet_started() -> None:
    ran: list[int] = []

    def piece(k: int) -> int:
        ran.append(k)
        time.sleep(0.05)
        return k

    def interrupt() -> None:
        raise _Interrupted

    with pytest.raises(_Interrupted):
        Workers(2).run([partial(piece, k) for k in range(8)], interrupt)
    time.sleep(0.2)
    assert len(ran) < 8


def _noisy(k: int) -> int:
    os.write(1, b"noise from a C level write\n")
    print("noise from print")
    return k * k


def test_what_a_piece_writes_to_standard_output_leaves_the_pipe_alone() -> None:
    workers = Workers(2)
    assert workers.apply(_noisy, (3,)) == 9
    assert workers.apply(_noisy, (4,)) == 16


def test_error_in_a_process_carries_its_traceback_as_a_note() -> None:
    with pytest.raises(ValueError, match="piece 1") as caught:
        Workers(2).solve(_fail_on_one, [(0,), (1,)])
    notes = getattr(caught.value, "__notes__", [])
    assert any("worker process" in note and "_fail_on_one" in note for note in notes)
