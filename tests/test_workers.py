import os
import threading
import time
from functools import partial

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
