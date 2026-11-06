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
