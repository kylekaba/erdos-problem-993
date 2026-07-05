#!/usr/bin/env python3
"""
Utilities for experimenting with Erdős problem 993 and the candidate
root/leaf-transfer lemmas discussed in the handoff memo.

The code is deliberately plain Python with exact integer/rational arithmetic.
NetworkX is optional and used only for non-isomorphic tree enumeration.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Iterable, List, Sequence, Tuple, Dict, Optional

Poly = List[int]
Edge = Tuple[int, int]


# ---------------------------------------------------------------------------
# Basic polynomial arithmetic over Z[x]
# ---------------------------------------------------------------------------

def trim(p: Poly) -> Poly:
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def coeff(p: Sequence[int], k: int) -> int:
    return p[k] if 0 <= k < len(p) else 0


def add(p: Sequence[int], q: Sequence[int]) -> Poly:
    n = max(len(p), len(q))
    return trim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def sub(p: Sequence[int], q: Sequence[int]) -> Poly:
    n = max(len(p), len(q))
    return trim([(p[i] if i < len(p) else 0) - (q[i] if i < len(q) else 0) for i in range(n)])


def mul(p: Sequence[int], q: Sequence[int]) -> Poly:
    if not p or not q:
        return [0]
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                if y:
                    r[i + j] += x * y
    return trim(r)


def shift(p: Sequence[int], s: int = 1) -> Poly:
    if p == [0] or not p:
        return [0]
    return [0] * s + list(p)


def pow_poly(p: Sequence[int], e: int) -> Poly:
    if e < 0:
        raise ValueError("negative polynomial exponent")
    base = list(p)
    out = [1]
    while e:
        if e & 1:
            out = mul(out, base)
        base = mul(base, base)
        e >>= 1
    return out


def prod(polys: Iterable[Sequence[int]]) -> Poly:
    out = [1]
    for p in polys:
        out = mul(out, p)
    return out


# ---------------------------------------------------------------------------
# Independence polynomials of forests by exact tree DP
# ---------------------------------------------------------------------------

def normalize_edges(n: int, edges: Iterable[Edge]) -> List[Edge]:
    out = []
    for u, v in edges:
        if not (0 <= u < n and 0 <= v < n):
            raise ValueError(f"edge {(u, v)} outside vertex range 0..{n-1}")
        if u == v:
            raise ValueError("loops are not allowed")
        out.append((u, v))
    return out


def adjacency(n: int, edges: Iterable[Edge]) -> List[List[int]]:
    adj = [[] for _ in range(n)]
    for u, v in normalize_edges(n, edges):
        adj[u].append(v)
        adj[v].append(u)
    return adj


def induced_edges(edges: Iterable[Edge], keep: Sequence[int]) -> Tuple[int, List[Edge], Dict[int, int]]:
    """Return relabeled induced subgraph on vertices in keep."""
    keep_set = set(keep)
    mp = {v: i for i, v in enumerate(sorted(keep_set))}
    new_edges = []
    for u, v in edges:
        if u in keep_set and v in keep_set:
            new_edges.append((mp[u], mp[v]))
    return len(mp), new_edges, mp


def delete_vertices(n: int, edges: Iterable[Edge], delete: Iterable[int]) -> Tuple[int, List[Edge], Dict[int, int]]:
    delete_set = set(delete)
    keep = [v for v in range(n) if v not in delete_set]
    return induced_edges(edges, keep)


def independence_poly_forest(n: int, edges: Iterable[Edge]) -> Poly:
    """Exact independence polynomial of a forest. Raises if the graph is not a forest."""
    edges = normalize_edges(n, edges)
    adj = adjacency(n, edges)
    seen = [False] * n
    whole = [1]

    edge_count = len(edges)
    # A forest on n vertices with c components has edge_count = n-c.
    # We also verify no cycles while traversing.
    for start in range(n):
        if seen[start]:
            continue
        parent = {start: -1}
        order = [start]
        seen[start] = True
        q = deque([start])
        while q:
            v = q.popleft()
            for u in adj[v]:
                if u == parent.get(v, -2):
                    continue
                if seen[u]:
                    raise ValueError("graph is not a forest: cycle detected")
                seen[u] = True
                parent[u] = v
                order.append(u)
                q.append(u)

        # DP for this rooted tree component.
        p0: Dict[int, Poly] = {v: [1] for v in order}      # v not selected
        p1: Dict[int, Poly] = {v: [0, 1] for v in order}   # v selected
        for v in reversed(order):
            for u in adj[v]:
                if parent.get(u) == v:
                    p1[v] = mul(p1[v], p0[u])
                    p0[v] = mul(p0[v], add(p0[u], p1[u]))
        comp_poly = add(p0[start], p1[start])
        whole = mul(whole, comp_poly)
    return trim(whole)


def independence_poly_after_deleting(n: int, edges: Iterable[Edge], delete: Iterable[int]) -> Poly:
    nn, ee, _ = delete_vertices(n, edges, delete)
    return independence_poly_forest(nn, ee)


def closed_neighborhood(n: int, edges: Iterable[Edge], root: int) -> List[int]:
    adj = adjacency(n, edges)
    return sorted(set([root] + adj[root]))


# ---------------------------------------------------------------------------
# Bush family and mixed bushes
# ---------------------------------------------------------------------------

def bush_edges(c: int, m: int) -> Tuple[int, List[Edge]]:
    """B(c,m): hub 0, m branch vertices, each with c pendant P2's."""
    edges: List[Edge] = []
    nid = 1
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


def bush_poly_closed(c: int, m: int) -> Poly:
    """I(B(c,m);x) = [(1+2x)^c + x(1+x)^c]^m + x(1+2x)^{mc}."""
    arm = add(pow_poly([1, 2], c), shift(pow_poly([1, 1], c)))
    return add(pow_poly(arm, m), shift(pow_poly([1, 2], m * c)))


def mixed_bush_edges(cs: Sequence[int]) -> Tuple[int, List[Edge]]:
    """Hub with arm i having c_i pendant P2's."""
    edges: List[Edge] = []
    nid = 1
    for c in cs:
        vi = nid
        nid += 1
        edges.append((0, vi))
        for _ in range(c):
            a, b = nid, nid + 1
            nid += 2
            edges.append((vi, a))
            edges.append((a, b))
    return nid, edges


def mixed_bush_poly_closed(cs: Sequence[int]) -> Poly:
    """Closed form for a hub with heterogeneous arm sizes c_i."""
    out_phase = [1]
    total_c = 0
    for c in cs:
        total_c += c
        arm = add(pow_poly([1, 2], c), shift(pow_poly([1, 1], c)))
        out_phase = mul(out_phase, arm)
    in_phase = shift(pow_poly([1, 2], total_c))
    return add(out_phase, in_phase)


# ---------------------------------------------------------------------------
# Sequence diagnostics
# ---------------------------------------------------------------------------

def mu_sequence(p: Sequence[int]) -> List[Fraction]:
    return [Fraction((k + 1) * p[k + 1], p[k]) for k in range(len(p) - 1) if p[k] != 0]


def max_mu_increment(p: Sequence[int]) -> Tuple[int, Fraction]:
    mu = mu_sequence(p)
    if len(mu) < 2:
        return -1, Fraction(0)
    k = max(range(len(mu) - 1), key=lambda j: mu[j + 1] - mu[j])
    return k, mu[k + 1] - mu[k]


def hinge_violations(p: Sequence[int]) -> List[int]:
    mu = mu_sequence(p)
    return [k for k in range(len(mu) - 1) if mu[k] <= k + 1 and mu[k + 1] > k + 2]


def is_unimodal(p: Sequence[int]) -> bool:
    seen_drop = False
    for i in range(len(p) - 1):
        if p[i + 1] < p[i]:
            seen_drop = True
        elif seen_drop and p[i + 1] > p[i]:
            return False
    return True


def strict_descent_indices(p: Sequence[int]) -> List[int]:
    return [k for k in range(len(p) - 1) if p[k + 1] < p[k]]


def first_differences(p: Sequence[int]) -> List[int]:
    return [coeff(p, k) - coeff(p, k - 1) for k in range(len(p) + 1)]


# ---------------------------------------------------------------------------
# RD / IRD conditions
# ---------------------------------------------------------------------------

def rd_failures(A: Sequence[int], B: Sequence[int]) -> List[Tuple[int, int, int, int]]:
    """Return failures of a_k<a_{k-1} => a_{k+1}+b_k <= a_k+b_{k-1}.

    Each tuple is (k, lhs, rhs, a_k_minus_a_kminus1).
    """
    maxk = max(len(A), len(B)) + 2
    bad = []
    for k in range(maxk):
        if coeff(A, k) < coeff(A, k - 1):
            lhs = coeff(A, k + 1) + coeff(B, k)
            rhs = coeff(A, k) + coeff(B, k - 1)
            if lhs > rhs:
                bad.append((k, lhs, rhs, coeff(A, k) - coeff(A, k - 1)))
    return bad


def has_rd(A: Sequence[int], B: Sequence[int]) -> bool:
    return not rd_failures(A, B)


def ird_failures(A: Sequence[int], B: Sequence[int]) -> List[Tuple[int, int, int, int]]:
    """IRD: for every interval [s,t], a_t<a_{s-1} => c_{t+1}<=c_s,
    where C=A+xB. Tuple: (s,t,lhs,rhs). Uses coefficient 0 outside support.
    """
    # Need only finite range around supports; after both A and B vanish, inequalities are trivial.
    maxidx = max(len(A), len(B)) + 2
    bad = []
    for s in range(1, maxidx + 1):
        for t in range(s, maxidx + 1):
            if coeff(A, t) < coeff(A, s - 1):
                lhs = coeff(A, t + 1) + coeff(B, t)
                rhs = coeff(A, s) + coeff(B, s - 1)
                if lhs > rhs:
                    bad.append((s, t, lhs, rhs))
    return bad


def has_ird(A: Sequence[int], B: Sequence[int]) -> bool:
    return not ird_failures(A, B)


def abstract_product_closure_counterexample() -> Dict[str, object]:
    """A small sequence-level counterexample to arbitrary RD product closure."""
    A = [1, 7, 8, 8, 8, 7]
    B = [1, 2, 4, 8, 8, 7]
    C = [1, 6, 8]
    D = [1, 4, 7]
    AC = mul(A, C)
    BD = mul(B, D)
    return {
        "A": A, "B": B, "C": C, "D": D,
        "AC": AC, "BD": BD,
        "rd_AB": has_rd(A, B),
        "rd_CD": has_rd(C, D),
        "rd_product": has_rd(AC, BD),
        "rd_product_failures": rd_failures(AC, BD),
    }


# ---------------------------------------------------------------------------
# Lemma tests on actual trees/forests
# ---------------------------------------------------------------------------

@dataclass
class LemmaFailure:
    lemma: str
    n: int
    edges: List[Edge]
    root_or_leaf: int
    k: int
    details: str


def root_pair_polys(n: int, edges: Sequence[Edge], root: int) -> Tuple[Poly, Poly, Poly]:
    """Return P=I(F), A=I(F-root), B=I(F-N[root])."""
    P = independence_poly_forest(n, edges)
    A = independence_poly_after_deleting(n, edges, [root])
    B = independence_poly_after_deleting(n, edges, closed_neighborhood(n, edges, root))
    return P, A, B


def check_compact_target(n: int, edges: Sequence[Edge]) -> List[LemmaFailure]:
    """Check I(F-r)_k<I(F-r)_{k-1} => I(F)_{k+1}<=I(F)_k for every root."""
    bad: List[LemmaFailure] = []
    for r in range(n):
        P, A, _B = root_pair_polys(n, edges, r)
        maxk = max(len(P), len(A)) + 1
        for k in range(maxk):
            if coeff(A, k) < coeff(A, k - 1) and coeff(P, k + 1) > coeff(P, k):
                bad.append(LemmaFailure(
                    "compact_root_deletion_dominance", n, list(edges), r, k,
                    f"A_k={coeff(A,k)} < A_k-1={coeff(A,k-1)} but P_k+1={coeff(P,k+1)} > P_k={coeff(P,k)}"
                ))
    return bad


def check_single_root_interlacing(n: int, edges: Sequence[Edge]) -> List[LemmaFailure]:
    """Check I(F)_{k+1}>I(F)_k => I(F-r)_k>=I(F-r)_{k-1}."""
    bad: List[LemmaFailure] = []
    for r in range(n):
        P, A, _B = root_pair_polys(n, edges, r)
        maxk = max(len(P), len(A)) + 1
        for k in range(maxk):
            if coeff(P, k + 1) > coeff(P, k) and coeff(A, k) < coeff(A, k - 1):
                bad.append(LemmaFailure(
                    "single_root_interlacing", n, list(edges), r, k,
                    f"P_k+1={coeff(P,k+1)} > P_k={coeff(P,k)} but A_k={coeff(A,k)} < A_k-1={coeff(A,k-1)}"
                ))
    return bad


def check_root_ird(n: int, edges: Sequence[Edge]) -> List[LemmaFailure]:
    """For every root r, test IRD for A=I(F-r), B=I(F-N[r])."""
    bad: List[LemmaFailure] = []
    for r in range(n):
        _P, A, B = root_pair_polys(n, edges, r)
        fails = ird_failures(A, B)
        for s, t, lhs, rhs in fails[:5]:
            bad.append(LemmaFailure(
                "root_pair_IRD", n, list(edges), r, t,
                f"interval s={s}, t={t}: lhs={lhs} > rhs={rhs}"
            ))
    return bad


def leaves(n: int, edges: Sequence[Edge]) -> List[int]:
    adj = adjacency(n, edges)
    return [v for v in range(n) if len(adj[v]) <= 1]  # isolated vertices count as leaves in forests


def check_leaf_no_hidden_strict_descent(n: int, edges: Sequence[Edge]) -> List[LemmaFailure]:
    """Check strict leaf transfer:
    If I(T)_{k+1}<I(T)_k, then both T-leaf and T-{leaf,neighbor} are not rising at k.
    For isolated vertex leaf, T-N[leaf]=T-leaf.
    """
    bad: List[LemmaFailure] = []
    P = independence_poly_forest(n, edges)
    adj = adjacency(n, edges)
    for ell in leaves(n, edges):
        if adj[ell]:
            u = adj[ell][0]
            del_leaf = [ell]
            del_pair = [ell, u]
        else:
            del_leaf = [ell]
            del_pair = [ell]
        A = independence_poly_after_deleting(n, edges, del_leaf)
        H = independence_poly_after_deleting(n, edges, del_pair)
        maxk = max(len(P), len(A), len(H)) + 1
        for k in range(maxk):
            if coeff(P, k + 1) < coeff(P, k):
                if coeff(A, k + 1) > coeff(A, k):
                    bad.append(LemmaFailure(
                        "leaf_delete_hidden_rise", n, list(edges), ell, k,
                        f"P descends {coeff(P,k+1)}<{coeff(P,k)} but T-leaf rises {coeff(A,k+1)}>{coeff(A,k)}"
                    ))
                if coeff(H, k + 1) > coeff(H, k):
                    bad.append(LemmaFailure(
                        "leaf_pair_delete_hidden_rise", n, list(edges), ell, k,
                        f"P descends {coeff(P,k+1)}<{coeff(P,k)} but T-{{leaf,neighbor}} rises {coeff(H,k+1)}>{coeff(H,k)}"
                    ))
    return bad


# ---------------------------------------------------------------------------
# Optional non-isomorphic tree enumeration via NetworkX
# ---------------------------------------------------------------------------

def nx_tree_edges(n: int):
    if n == 1:
        yield []
        return
    try:
        import networkx as nx
    except Exception as e:  # pragma: no cover
        raise RuntimeError("NetworkX is required for non-isomorphic tree enumeration") from e
    for T in nx.generators.nonisomorphic_trees(n):
        yield sorted((min(u, v), max(u, v)) for u, v in T.edges())


def test_all_nonisomorphic_trees(nmax: int = 12, tests: Sequence[str] = ("compact", "interlacing", "ird", "leaf")) -> List[LemmaFailure]:
    bad: List[LemmaFailure] = []
    for n in range(1, nmax + 1):
        for edges in nx_tree_edges(n):
            if "compact" in tests:
                bad.extend(check_compact_target(n, edges))
            if "interlacing" in tests:
                bad.extend(check_single_root_interlacing(n, edges))
            if "ird" in tests:
                bad.extend(check_root_ird(n, edges))
            if "leaf" in tests:
                bad.extend(check_leaf_no_hidden_strict_descent(n, edges))
            if bad:
                return bad
    return bad


# ---------------------------------------------------------------------------
# Pretty reports
# ---------------------------------------------------------------------------

def summarize_poly(p: Sequence[int], max_terms: int = 12) -> str:
    if len(p) <= max_terms:
        return str(list(p))
    return f"{list(p[:max_terms])} ... {list(p[-3:])} (len={len(p)})"


def report_bush(c: int, m: int) -> Dict[str, object]:
    p = bush_poly_closed(c, m)
    k, inc = max_mu_increment(p)
    mu = mu_sequence(p)
    return {
        "c": c,
        "m": m,
        "n": m * (2 * c + 1) + 1,
        "alpha": len(p) - 1,
        "max_increment_k": k,
        "mu_k": mu[k] if k >= 0 else None,
        "mu_k_plus_1": mu[k + 1] if k >= 0 and k + 1 < len(mu) else None,
        "increment": inc,
        "increment_float": float(inc),
        "exceeds_1": inc > 1,
        "unimodal": is_unimodal(p),
        "hinge_violations": hinge_violations(p),
    }


if __name__ == "__main__":
    # Minimal smoke test when run directly.
    for c, m in [(8, 14), (8, 15), (10, 36), (12, 80)]:
        r = report_bush(c, m)
        print(f"B({c},{m}) n={r['n']} alpha={r['alpha']} k={r['max_increment_k']} inc={r['increment_float']:.6f} unimodal={r['unimodal']} hinge_bad={r['hinge_violations']}")
