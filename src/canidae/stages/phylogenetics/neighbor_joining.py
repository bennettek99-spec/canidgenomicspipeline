"""Neighbor-joining tree construction with bootstrap support (Saitou & Nei 1987).

A dependency-light NJ implementation over a distance matrix, producing a Newick string with
branch lengths. Bootstrap support is obtained by resampling SNP sites with replacement,
rebuilding the tree, and counting how often each bipartition of the original tree recurs.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Clade:
    """A node in the tree. Leaves have a ``name`` and no children."""

    name: str | None
    length: float
    tips: frozenset[str]
    children: list[Clade] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children


def neighbor_joining(dist: np.ndarray, labels: list[str]) -> Clade:
    """Build an NJ tree from a square distance matrix; return the (unrooted) root clade."""
    n = len(labels)
    if n < 2:
        raise ValueError("neighbor_joining needs >= 2 taxa")
    clades: list[Clade] = [Clade(name=lab, length=0.0, tips=frozenset([lab])) for lab in labels]
    d: dict[int, dict[int, float]] = {i: {} for i in range(n)}
    for i in range(n):
        for j in range(n):
            if i != j:
                d[i][j] = float(dist[i][j])
    active = list(range(n))

    while len(active) > 2:
        m = len(active)
        r = {i: sum(d[i][j] for j in active if j != i) for i in active}
        best_pair, best_q = None, np.inf
        for ai in range(m):
            for bi in range(ai + 1, m):
                a, b = active[ai], active[bi]
                q = (m - 2) * d[a][b] - r[a] - r[b]
                if q < best_q:
                    best_q, best_pair = q, (a, b)
        assert best_pair is not None
        a, b = best_pair
        d_ab = d[a][b]
        da = 0.5 * d_ab + (r[a] - r[b]) / (2 * (m - 2))
        db = d_ab - da
        clades[a].length = da
        clades[b].length = db

        u = Clade(
            name=None,
            length=0.0,
            tips=clades[a].tips | clades[b].tips,
            children=[clades[a], clades[b]],
        )
        uid = len(clades)
        clades.append(u)
        d[uid] = {}
        for k in active:
            if k in (a, b):
                continue
            duk = 0.5 * (d[a][k] + d[b][k] - d_ab)
            d[uid][k] = duk
            d[k][uid] = duk
        active.remove(a)
        active.remove(b)
        active.append(uid)

    a, b = active
    d_ab = d[a][b]
    clades[a].length = clades[b].length = d_ab / 2
    return Clade(
        name=None, length=0.0, tips=clades[a].tips | clades[b].tips, children=[clades[a], clades[b]]
    )


# --------------------------------------------------------------------------------------
# Bipartitions & support
# --------------------------------------------------------------------------------------


def _iter_internal(clade: Clade):
    if not clade.is_leaf:
        yield clade
        for child in clade.children:
            yield from _iter_internal(child)


def bipartitions(root: Clade) -> set[frozenset[str]]:
    """Canonical non-trivial bipartitions of a tree (root-independent)."""
    all_tips = root.tips
    ref = min(all_tips)  # fixed reference for canonicalization
    splits: set[frozenset[str]] = set()
    for node in _iter_internal(root):
        side = node.tips if ref not in node.tips else (all_tips - node.tips)
        if 2 <= len(side) <= len(all_tips) - 1:
            splits.add(frozenset(side))
    return splits


def bootstrap_support(
    gn: np.ndarray,
    labels: list[str],
    distance_fn,
    *,
    n_boot: int = 100,
    seed: int = 0,
    blocks: np.ndarray | None = None,
) -> tuple[Clade, dict[frozenset[str], float]]:
    """Return the NJ tree on the full data plus per-bipartition bootstrap support in [0,1].

    ``gn`` is a (n_sites, n_taxa) matrix; ``distance_fn(gn)`` returns a (n_taxa, n_taxa)
    distance matrix. Sites are resampled with replacement unless block labels are supplied,
    in which case complete linkage blocks are resampled together.
    """
    main = neighbor_joining(distance_fn(gn), labels)
    target = bipartitions(main)
    counts: dict[frozenset[str], int] = dict.fromkeys(target, 0)

    n_sites = gn.shape[0]
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        if blocks is None:
            cols = rng.integers(0, n_sites, n_sites)
        else:
            block_array = np.asarray(blocks)
            unique = np.unique(block_array)
            sampled = rng.choice(unique, size=unique.size, replace=True)
            cols = np.concatenate([np.flatnonzero(block_array == value) for value in sampled])
        boot = bipartitions(neighbor_joining(distance_fn(gn[cols]), labels))
        for split in target:
            if split in boot:
                counts[split] += 1
    support = {s: counts[s] / n_boot for s in target} if n_boot else {}
    return main, support


# --------------------------------------------------------------------------------------
# Newick
# --------------------------------------------------------------------------------------


def to_newick(root: Clade, support: dict[frozenset[str], float] | None = None) -> str:
    """Render a Newick string; internal nodes labeled with bootstrap support (%) if given."""
    all_tips = root.tips
    ref = min(all_tips)
    support = support or {}

    def label_for(node: Clade) -> str:
        side = node.tips if ref not in node.tips else (all_tips - node.tips)
        val = support.get(frozenset(side))
        return "" if val is None else f"{val * 100:.0f}"

    def render(node: Clade, *, is_root: bool = False) -> str:
        if node.is_leaf:
            return f"{node.name}:{max(node.length, 0.0):.6f}"
        inner = ",".join(render(c) for c in node.children)
        if is_root:
            return f"({inner});"
        return f"({inner}){label_for(node)}:{max(node.length, 0.0):.6f}"

    return render(root, is_root=True)
