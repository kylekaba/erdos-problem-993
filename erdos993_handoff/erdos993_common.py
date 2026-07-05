#!/usr/bin/env python3
"""Shared exact-arithmetic utilities for Erdős problem 993 experiments.

All polynomials are coefficient lists [a_0, a_1, ...] with Python ints.
Graphs are forests represented by (n, edges), with vertices 0..n-1.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb
from typing import Iterable, List, Sequence, Tuple, Set

Poly = List[int]
Edge = Tuple[int, int]


def trim(p: Sequence[int]) -> Poly:
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p or [0]


def coeff(p: Sequence[int], k: int) -> int:
    return p[k] if 0 <= k < len(p) else 0


def poly_add(a: Sequence[int], b: Sequence[int]) -> Poly:
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])


def poly_sub(a: Sequence[int], b: Sequence[int]) -> Poly:
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(n)])


def poly_shift(a: Sequence[int], s: int = 1) -> Poly:
    if not a or a == [0]:
        return [0]
    return [0] * s + list(a)


def poly_mul(a: Sequence[int], b: Sequence[int]) -> Poly:
    if not a or not b or a == [0] or b == [0]:
        return [0]
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] += x * y
    return trim(r)


def poly_pow(p: Sequence[int], e: int) -> Poly:
    r = [1]
    p = list(p)
    while e:
        if e & 1:
            r = poly_mul(r, p)
        p = poly_mul(p, p)
        e >>= 1
    return r


def product(polys: Iterable[Sequence[int]]) -> Poly:
    r = [1]
    for p in polys:
        r = poly_mul(r, p)
    return r


def adjacency(n: int, edges: Iterable[Edge]) -> List[List[int]]:
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def induced_subforest(n: int, edges: Iterable[Edge], remove: Set[int]) -> Tuple[int, List[Edge], dict[int, int]]:
    """Return induced forest after deleting remove, with vertices relabeled 0..n'-1."""
    keep = [v for v in range(n) if v not in remove]
    relabel = {v: i for i, v in enumerate(keep)}
    new_edges: List[Edge] = []
    for u, v in edges:
        if u not in remove and v not in remove:
            new_edges.append((relabel[u], relabel[v]))
    return len(keep), new_edges, relabel


def closed_neighborhood(n: int, edges: Iterable[Edge], r: int) -> Set[int]:
    adj = adjacency(n, edges)
    return {r, *adj[r]}


def leaves(n: int, edges: Iterable[Edge]) -> List[int]:
    if n == 1:
        return [0]
    adj = adjacency(n, edges)
    return [v for v in range(n) if len(adj[v]) == 1]


def independence_polynomial(n: int, edges: Iterable[Edge]) -> Poly:
    """Exact independence polynomial of a forest via include/exclude tree DP."""
    edges = list(edges)
    adj = adjacency(n, edges)
    seen = [False] * n
    total = [1]

    for root in range(n):
        if seen[root]:
            continue
        parent = {root: -1}
        order = [root]
        seen[root] = True
        i = 0
        while i < len(order):
            v = order[i]
            i += 1
            for u in adj[v]:
                if not seen[u]:
                    seen[u] = True
                    parent[u] = v
                    order.append(u)

        # p0[v]: polynomial for subtree at v excluding v.
        # p1[v]: polynomial for subtree at v including v.
        p0 = {v: [1] for v in order}
        p1 = {v: [0, 1] for v in order}
        for v in reversed(order):
            for u in adj[v]:
                if parent.get(u) == v:
                    p1[v] = poly_mul(p1[v], p0[u])
                    p0[v] = poly_mul(p0[v], poly_add(p0[u], p1[u]))
        comp = poly_add(p0[root], p1[root])
        total = poly_mul(total, comp)
    return trim(total)


def path_polynomial(n_vertices: int) -> Poly:
    """I(P_n;x) for a path on n vertices."""
    return [comb(n_vertices - k + 1, k) for k in range(n_vertices // 2 + 1)]


def bush_tree(c: int, m: int) -> Tuple[int, List[Edge]]:
    """B(c,m): hub; m arms; each arm carries c pendant paths of length 2."""
    edges: List[Edge] = []
    nid = 1  # hub is 0
    for _ in range(m):
        vi = nid
        nid += 1
        edges.append((0, vi))
        for _ in range(c):
            a, b = nid, nid + 1
            nid += 2
            edges.append((vi, a))
            edges.append((a, b))
    return nid, edges


def bush_closed_polynomial(c: int, m: int) -> Poly:
    """Closed form: [(1+2x)^c + x(1+x)^c]^m + x(1+2x)^(mc)."""
    arm = poly_add(poly_pow([1, 2], c), poly_shift(poly_pow([1, 1], c)))
    return poly_add(poly_pow(arm, m), poly_shift(poly_pow([1, 2], m * c)))


def mixed_bush_closed_polynomial(cs: Sequence[int]) -> Poly:
    """Heterogeneous bush: one hub; arm i has c_i pendant length-2 paths."""
    arms = []
    total_c = 0
    for c in cs:
        total_c += c
        arms.append(poly_add(poly_pow([1, 2], c), poly_shift(poly_pow([1, 1], c))))
    return poly_add(product(arms), poly_shift(poly_pow([1, 2], total_c)))


def mu_sequence(poly: Sequence[int]) -> List[Fraction]:
    return [Fraction((k + 1) * poly[k + 1], poly[k]) for k in range(len(poly) - 1) if poly[k] != 0]


def max_mu_increment(poly: Sequence[int]) -> Tuple[int, Fraction]:
    mu = mu_sequence(poly)
    best_k = max(range(len(mu) - 1), key=lambda k: mu[k + 1] - mu[k])
    return best_k, mu[best_k + 1] - mu[best_k]


def is_unimodal(poly: Sequence[int]) -> bool:
    down = False
    for i in range(len(poly) - 1):
        if poly[i + 1] < poly[i]:
            down = True
        elif down and poly[i + 1] > poly[i]:
            return False
    return True


def hinge_violations(poly: Sequence[int]) -> List[int]:
    """Return k with mu_k <= k+1 but mu_{k+1} > k+2."""
    mu = mu_sequence(poly)
    return [k for k in range(len(mu) - 1) if mu[k] <= k + 1 and mu[k + 1] > k + 2]


def is_log_concave(seq: Sequence[int]) -> bool:
    return all(seq[i] * seq[i] >= seq[i - 1] * seq[i + 1] for i in range(1, len(seq) - 1))
