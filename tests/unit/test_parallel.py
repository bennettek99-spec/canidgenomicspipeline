from __future__ import annotations

import pytest

from canidae.core.parallel import map_parallel, scatter_intervals


def _square(x: int) -> int:
    return x * x


def _boom(_: int) -> int:
    raise RuntimeError("task failed")


def test_map_parallel_empty_input() -> None:
    assert map_parallel(_square, [], max_workers=4) == []


def test_map_parallel_single_item_skips_pool() -> None:
    assert map_parallel(_square, [3], max_workers=8) == [9]


def test_map_parallel_workers_clamped_to_item_count() -> None:
    # Asking for 64 workers over 3 items must still produce every result.
    assert map_parallel(_square, [2, 3, 4], max_workers=64) == [4, 9, 16]


def test_map_parallel_preserves_order_by_default() -> None:
    items = list(range(50))
    results = map_parallel(lambda x: x * 10, items, max_workers=8)
    assert results == [x * 10 for x in items]


def test_map_parallel_unordered_returns_same_multiset() -> None:
    items = list(range(32))
    results = map_parallel(_square, items, max_workers=6, ordered=False)
    assert sorted(results) == sorted(x * x for x in items)


def test_map_parallel_exception_propagates_from_first_failure() -> None:
    with pytest.raises(RuntimeError, match="task failed"):
        map_parallel(lambda x: _boom(x), [1], max_workers=2)


def test_map_parallel_process_backend() -> None:
    assert map_parallel(_square, [5, 6, 7], backend="process", max_workers=2) == [25, 36, 49]


def test_scatter_intervals_splits_fixed_windows() -> None:
    intervals = scatter_intervals([("chr1", 25_000_000)], window=10_000_000)
    assert [(iv.contig, iv.start, iv.end) for iv in intervals] == [
        ("chr1", 0, 10_000_000),
        ("chr1", 10_000_000, 20_000_000),
        ("chr1", 20_000_000, 25_000_000),
    ]


def test_scatter_intervals_exact_multiple_has_no_tail() -> None:
    intervals = scatter_intervals([("chr2", 20)], window=10)
    assert [(iv.start, iv.end) for iv in intervals] == [(0, 10), (10, 20)]


def test_scatter_intervals_multiple_contigs_and_empty_input() -> None:
    intervals = scatter_intervals([("chrA", 15), ("chrB", 5)], window=10)
    assert [(iv.contig, iv.start, iv.end) for iv in intervals] == [
        ("chrA", 0, 10),
        ("chrA", 10, 15),
        ("chrB", 0, 5),
    ]
    assert scatter_intervals([], window=10) == []
