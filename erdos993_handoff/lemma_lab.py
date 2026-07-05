#!/usr/bin/env python3
"""Search/verify proposed local lemmas for Erdős 993 on exact tree data.

Default mode checks all non-isomorphic trees up to --max-n using networkx.
The script prints the first counterexample if found, otherwise a summary.

Checks implemented:
  RD        : I(F-r)_k < I(F-r)_{k-1} => I(F)_{k+1} <= I(F)_k.
  INTERLACE : I(F)_{k+1} > I(F)_k => I(F-r)_k >= I(F-r)_{k-1}.
  IRD       : interval version for pair (F-r, F-N[r]).
  LEAF      : strict descent in T cannot be hidden by deleting a leaf or leaf+neighbor.
  TWO_NHR   : h_{k+1}>h_k>=h_{k-1} => h_{k+1}+j_k > h_{k-1}+j_{k-1}
              for pairs (H,J)=(F-r,F-N[r]).

These are verification tools, not proofs.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Iterable, List, Tuple

import networkx as nx

from erdos993_common import (
    adjacency,
    closed_neighborhood,
    coeff,
    induced_subforest,
    independence_polynomial,
    leaves,
)


@dataclass
class Failure:
    check: str
    n: int
    edges: List[Tuple[int, int]]
    root_or_leaf: int | Tuple[int, int]
    k_or_interval: int | Tuple[int, int]
    detail: str


def poly_after_delete(n: int, edges: List[Tuple[int, int]], remove: set[int]) -> List[int]:
    n2, e2, _ = induced_subforest(n, edges, remove)
    return independence_polynomial(n2, e2)


def check_rd(n: int, edges: List[Tuple[int, int]], r: int) -> Failure | None:
    p = independence_polynomial(n, edges)
    a = poly_after_delete(n, edges, {r})
    b = poly_after_delete(n, edges, closed_neighborhood(n, edges, r))
    maxk = max(len(a), len(b), len(p)) + 1
    for k in range(1, maxk):
        if coeff(a, k) < coeff(a, k - 1) and coeff(p, k + 1) > coeff(p, k):
            return Failure("RD", n, edges, r, k, f"a_k={coeff(a,k)} < a_k-1={coeff(a,k-1)} but p_k+1={coeff(p,k+1)} > p_k={coeff(p,k)}; b_k={coeff(b,k)}, b_k-1={coeff(b,k-1)}")
    return None


def check_interlace(n: int, edges: List[Tuple[int, int]], r: int) -> Failure | None:
    p = independence_polynomial(n, edges)
    a = poly_after_delete(n, edges, {r})
    maxk = max(len(a), len(p)) + 1
    for k in range(1, maxk):
        if coeff(p, k + 1) > coeff(p, k) and coeff(a, k) < coeff(a, k - 1):
            return Failure("INTERLACE", n, edges, r, k, f"p_k+1={coeff(p,k+1)} > p_k={coeff(p,k)} but a_k={coeff(a,k)} < a_k-1={coeff(a,k-1)}")
    return None


def check_ird(n: int, edges: List[Tuple[int, int]], r: int) -> Failure | None:
    # IRD for A=F-r, B=F-N[r], C=A+xB=F.
    a = poly_after_delete(n, edges, {r})
    p = independence_polynomial(n, edges)
    maxdeg = max(len(a), len(p)) + 1
    for s in range(1, maxdeg):
        for t in range(s, maxdeg):
            if coeff(a, t) < coeff(a, s - 1) and coeff(p, t + 1) > coeff(p, s):
                return Failure("IRD", n, edges, r, (s, t), f"a_t={coeff(a,t)} < a_s-1={coeff(a,s-1)} but p_t+1={coeff(p,t+1)} > p_s={coeff(p,s)}")
    return None


def check_leaf(n: int, edges: List[Tuple[int, int]], leaf: int) -> Failure | None:
    adj = adjacency(n, edges)
    if n == 1:
        return None
    u = adj[leaf][0]
    p = independence_polynomial(n, edges)
    a = poly_after_delete(n, edges, {leaf})
    h = poly_after_delete(n, edges, {leaf, u})
    maxk = max(len(p), len(a), len(h)) + 1
    for k in range(0, maxk):
        if coeff(p, k + 1) < coeff(p, k):
            if coeff(a, k + 1) > coeff(a, k):
                return Failure("LEAF_DELETE", n, edges, (leaf, u), k, f"p descends: {coeff(p,k+1)} < {coeff(p,k)} but T-leaf rises: {coeff(a,k+1)} > {coeff(a,k)}")
            if coeff(h, k + 1) > coeff(h, k):
                return Failure("LINK_TRANSFER", n, edges, (leaf, u), k, f"p descends: {coeff(p,k+1)} < {coeff(p,k)} but T-leaf-neighbor rises: {coeff(h,k+1)} > {coeff(h,k)}")
    return None


def check_two_nhr(n: int, edges: List[Tuple[int, int]], r: int) -> Failure | None:
    h = poly_after_delete(n, edges, {r})
    j = poly_after_delete(n, edges, closed_neighborhood(n, edges, r))
    maxk = max(len(h), len(j)) + 1
    for k in range(1, maxk):
        if coeff(h, k + 1) > coeff(h, k) >= coeff(h, k - 1):
            lhs = coeff(h, k + 1) + coeff(j, k)
            rhs = coeff(h, k - 1) + coeff(j, k - 1)
            if lhs <= rhs:
                return Failure("TWO_NHR", n, edges, r, k, f"h_k+1={coeff(h,k+1)} > h_k={coeff(h,k)} >= h_k-1={coeff(h,k-1)}, but h_k+1+j_k={lhs} <= h_k-1+j_k-1={rhs}")
    return None


def iter_trees(max_n: int):
    for n in range(1, max_n + 1):
        if n == 1:
            yield n, []
            continue
        for T in nx.nonisomorphic_trees(n):
            edges = sorted((min(u, v), max(u, v)) for u, v in T.edges())
            yield n, edges


def run_checks(max_n: int, checks: set[str]) -> None:
    tree_count = 0
    rooted_count = 0
    leaf_count = 0
    for n, edges in iter_trees(max_n):
        tree_count += 1
        if "LEAF" in checks:
            for leaf in leaves(n, edges):
                leaf_count += 1
                fail = check_leaf(n, edges, leaf)
                if fail:
                    print_failure(fail)
                    return
        for r in range(n):
            rooted_count += 1
            for name, fn in [("RD", check_rd), ("INTERLACE", check_interlace), ("IRD", check_ird), ("TWO_NHR", check_two_nhr)]:
                if name in checks:
                    fail = fn(n, edges, r)
                    if fail:
                        print_failure(fail)
                        return
    print(f"OK: checked {tree_count} non-isomorphic trees up to n={max_n}; {rooted_count} rooted cases; {leaf_count} leaf cases; checks={sorted(checks)}")


def print_failure(f: Failure) -> None:
    print("FAIL")
    print(f"check: {f.check}")
    print(f"n: {f.n}")
    print(f"edges: {f.edges}")
    print(f"root_or_leaf: {f.root_or_leaf}")
    print(f"k_or_interval: {f.k_or_interval}")
    print(f"detail: {f.detail}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=12, help="largest tree size to enumerate")
    ap.add_argument("--checks", nargs="+", default=["RD", "INTERLACE", "IRD", "LEAF", "TWO_NHR"], choices=["RD", "INTERLACE", "IRD", "LEAF", "TWO_NHR"])
    args = ap.parse_args()
    run_checks(args.max_n, set(args.checks))


if __name__ == "__main__":
    main()
