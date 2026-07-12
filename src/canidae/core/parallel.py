"""Parallel-execution helpers for stage internals.

The executor parallelizes *across* stages. These helpers parallelize *within* a stage —
the common bioinformatics pattern of scattering work across genomic intervals or samples,
then gathering. External tools do the heavy lifting in subprocesses, so a thread pool is
the right default (I/O-bound orchestration, no GIL contention); switch to processes only
for pure-Python compute.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import Literal, TypeVar

from canidae.core.model import GenomicInterval

T = TypeVar("T")
R = TypeVar("R")

Backend = Literal["thread", "process"]


def map_parallel(
    func: Callable[[T], R],
    items: Iterable[T],
    *,
    max_workers: int = 4,
    backend: Backend = "thread",
    ordered: bool = True,
) -> list[R]:
    """Apply ``func`` to each item concurrently and collect the results.

    ``ordered`` preserves input order in the output (the default, so gather is
    deterministic). Exceptions propagate from the first failing task.
    """
    items = list(items)
    if not items:
        return []
    workers = max(1, min(max_workers, len(items)))
    if workers == 1:
        return [func(x) for x in items]

    pool_cls = ThreadPoolExecutor if backend == "thread" else ProcessPoolExecutor
    with pool_cls(max_workers=workers) as pool:
        if ordered:
            return list(pool.map(func, items))
        futures = [pool.submit(func, x) for x in items]
        return [f.result() for f in futures]


def scatter_intervals(
    contigs: Sequence[tuple[str, int]],
    *,
    window: int = 10_000_000,
) -> list[GenomicInterval]:
    """Split ``(contig, length)`` pairs into fixed-size scatter intervals.

    A window of ~10 Mb is a sensible default for per-region variant calling: small enough
    to parallelize widely, large enough to amortize tool start-up.
    """
    intervals: list[GenomicInterval] = []
    for contig, length in contigs:
        start = 0
        while start < length:
            end = min(start + window, length)
            intervals.append(GenomicInterval(contig=contig, start=start, end=end))
            start = end
    return intervals
