"""Disjoint sets of integers."""

from __future__ import annotations


class UnionFind:
    """Disjoint sets of the integers below ``count``, joined by ``union``; ``find`` names
    each set by one member, its root."""

    def __init__(self, count: int) -> None:
        self.parent = list(range(count))

    def find(self, i: int) -> int:
        """Root of the set that holds ``i``."""
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def union(self, a: int, b: int) -> None:
        """Join the set of ``b`` to that of ``a``, whose root stays."""
        self.parent[self.find(b)] = self.find(a)
