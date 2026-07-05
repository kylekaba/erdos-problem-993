#!/usr/bin/env python3
"""Verify the uploaded bush counterexamples in exact arithmetic.

Usage:
  python verify_bush_counterexamples.py --json /mnt/data/erdos993_mu_increment_counterexamples.json

Checks:
  * closed-form polynomial for every homogeneous B(c,m) entry in the JSON;
  * exact mu_k values and maximum mu increment;
  * unimodality and hinge implication for each witness;
  * independent tree-DP equality for small homogeneous witnesses with explicit edges.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from erdos993_common import (
    bush_closed_polynomial,
    bush_tree,
    hinge_violations,
    independence_polynomial,
    is_unimodal,
    max_mu_increment,
    mu_sequence,
)


def parse_fraction(s: str) -> Fraction:
    return Fraction(s)


def check_entry(name: str, rec: dict, do_dp_limit: int) -> None:
    c = rec.get("c")
    m = rec.get("m")
    if not isinstance(c, int) or not isinstance(m, int):
        print(f"SKIP {name}: no homogeneous (c,m) metadata")
        return

    poly = bush_closed_polynomial(c, m)
    alpha = len(poly) - 1
    n = m * (2 * c + 1) + 1
    assert rec["n"] == n, (name, "n mismatch", rec["n"], n)
    assert rec["alpha"] == alpha, (name, "alpha mismatch", rec["alpha"], alpha)

    if "i_sequence" in rec:
        seq_json = [int(x) for x in rec["i_sequence"]]
        assert seq_json == poly, (name, "i_sequence differs from closed form")

    mu = mu_sequence(poly)
    k = rec["violating_k"]
    assert parse_fraction(rec["mu_k"]) == mu[k], (name, "mu_k mismatch")
    assert parse_fraction(rec["mu_k_plus_1"]) == mu[k + 1], (name, "mu_k_plus_1 mismatch")

    best_k, best_inc = max_mu_increment(poly)
    inc = mu[k + 1] - mu[k]
    assert inc > 1, (name, "claimed increment is not > 1")
    assert best_k == k, (name, "violating_k is not max increment", k, best_k)
    assert rec.get("exceeds_1") is True
    assert is_unimodal(poly) == rec.get("unimodal"), (name, "unimodality flag mismatch")
    assert hinge_violations(poly) == rec.get("hinge_violations", []), (name, "hinge list mismatch")

    if n <= do_dp_limit:
        n2, edges = bush_tree(c, m)
        assert n2 == n
        dp_poly = independence_polynomial(n2, edges)
        assert dp_poly == poly, (name, "tree DP differs from closed form")
        dp_msg = "; DP=closed form"
    else:
        dp_msg = ""

    print(
        f"OK {name}: B({c},{m}), n={n}, alpha={alpha}, "
        f"k={k}, mu_inc={float(inc):.12f}, unimodal={is_unimodal(poly)}, "
        f"hinge_violations={hinge_violations(poly)}{dp_msg}"
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="/mnt/data/erdos993_mu_increment_counterexamples.json")
    ap.add_argument("--dp-limit", type=int, default=300, help="run full tree DP for homogeneous bushes with n <= this")
    args = ap.parse_args()

    data = json.loads(Path(args.json).read_text())
    for name, rec in data.items():
        check_entry(name, rec, args.dp_limit)


if __name__ == "__main__":
    main()
