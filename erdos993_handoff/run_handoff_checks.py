#!/usr/bin/env python3
"""
Run exact verification checks for the Erdős 993 handoff.

Examples:
  python run_handoff_checks.py --bush
  python run_handoff_checks.py --json ../erdos993_mu_increment_counterexamples.json
  python run_handoff_checks.py --small-trees 12
  python run_handoff_checks.py --abstract-product-counterexample
  python run_handoff_checks.py --all
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from fractions import Fraction

from erdos993_tools import (
    abstract_product_closure_counterexample,
    bush_poly_closed,
    coeff,
    has_ird,
    hinge_violations,
    independence_poly_forest,
    is_unimodal,
    max_mu_increment,
    mixed_bush_poly_closed,
    mu_sequence,
    report_bush,
    test_all_nonisomorphic_trees,
)


def fmt_frac(q):
    if q is None:
        return "None"
    return f"{q.numerator}/{q.denominator}" if q.denominator != 1 else str(q.numerator)


def run_bush():
    print("\n=== Homogeneous bush witnesses ===")
    for c, m in [(8, 14), (8, 15), (10, 36), (12, 80)]:
        r = report_bush(c, m)
        print(f"B({c},{m}): n={r['n']}, alpha={r['alpha']}, k={r['max_increment_k']}")
        print(f"  mu_k       = {float(r['mu_k']):.12f} = {fmt_frac(r['mu_k'])}")
        print(f"  mu_k+1     = {float(r['mu_k_plus_1']):.12f} = {fmt_frac(r['mu_k_plus_1'])}")
        print(f"  increment  = {r['increment_float']:.12f} = {fmt_frac(r['increment'])}")
        print(f"  inc > 1?   = {r['exceeds_1']}")
        print(f"  unimodal?  = {r['unimodal']}")
        print(f"  hinge violations = {r['hinge_violations']}")
        assert r["exceeds_1"]
        assert r["unimodal"]
        assert not r["hinge_violations"]

    print("\n=== Mixed smaller witness: eight arms c=7 and five arms c=8 ===")
    cs = [7] * 8 + [8] * 5
    p = mixed_bush_poly_closed(cs)
    k, inc = max_mu_increment(p)
    mu = mu_sequence(p)
    print(f"mixed bush: n={1 + sum(2*c+1 for c in cs)}, alpha={len(p)-1}, k={k}")
    print(f"  mu_k      = {float(mu[k]):.12f}")
    print(f"  mu_k+1    = {float(mu[k+1]):.12f}")
    print(f"  increment = {float(inc):.12f} = {fmt_frac(inc)}")
    print(f"  inc > 1?  = {inc > 1}")
    print(f"  unimodal? = {is_unimodal(p)}")
    print(f"  hinge violations = {hinge_violations(p)}")
    assert inc > 1
    assert is_unimodal(p)
    assert not hinge_violations(p)


def run_json(path):
    print(f"\n=== Verify JSON records against closed forms: {path} ===")
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    for name, rec in data.items():
        c = int(rec["c"])
        m = int(rec["m"])
        p = bush_poly_closed(c, m)
        json_seq = [int(x) for x in rec["i_sequence"]]
        print(f"{name}: B({c},{m}), n={rec['n']}, alpha={rec['alpha']}")
        assert p == json_seq, f"closed form differs from JSON i_sequence for {name}"
        k, inc = max_mu_increment(p)
        mu = mu_sequence(p)
        print(f"  max increment k={k}, increment={float(inc):.12f}")
        print(f"  JSON violating_k={rec['violating_k']}, JSON increment={rec['increment_float']}")
        assert k == int(rec["violating_k"])
        assert inc > 1
        assert is_unimodal(p) == bool(rec["unimodal"])
        assert hinge_violations(p) == list(rec["hinge_violations"])
        assert str(mu[k]) == rec["mu_k"]
        assert str(mu[k + 1]) == rec["mu_k_plus_1"]
    print("JSON verification passed.")


def run_abstract_counterexample():
    print("\n=== Abstract product-closure counterexample ===")
    ce = abstract_product_closure_counterexample()
    for key in ["A", "B", "C", "D", "AC", "BD"]:
        print(f"{key} = {ce[key]}")
    print(f"RD(A,B)       = {ce['rd_AB']}")
    print(f"RD(C,D)       = {ce['rd_CD']}")
    print(f"RD(AC,BD)     = {ce['rd_product']}")
    print(f"failures      = {ce['rd_product_failures']}")
    assert ce["rd_AB"] and ce["rd_CD"] and not ce["rd_product"]


def run_small_trees(nmax):
    print(f"\n=== Exhaustive non-isomorphic tree tests through n={nmax} ===")
    print("Tests: compact target, single-root interlacing, root-pair IRD, strict leaf no-hidden-rise")
    failures = test_all_nonisomorphic_trees(nmax=nmax)
    if not failures:
        print("No failures found.")
        return
    print(f"FOUND {len(failures)} failure(s); first failure:")
    f = failures[0]
    print(f)
    raise SystemExit(1)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bush", action="store_true", help="verify homogeneous and mixed bush witnesses")
    ap.add_argument("--json", help="verify uploaded JSON records against closed forms")
    ap.add_argument("--abstract-product-counterexample", action="store_true", help="show RD is not closed under arbitrary products")
    ap.add_argument("--small-trees", type=int, metavar="N", help="test proposed lemmas on all nonisomorphic trees up to N")
    ap.add_argument("--all", action="store_true", help="run a useful default suite")
    args = ap.parse_args(argv)

    if args.all:
        run_bush()
        default_json = os.path.join(os.path.dirname(__file__), "..", "erdos993_mu_increment_counterexamples.json")
        if os.path.exists(default_json):
            run_json(default_json)
        run_abstract_counterexample()
        run_small_trees(12)
        return

    did = False
    if args.bush:
        run_bush(); did = True
    if args.json:
        run_json(args.json); did = True
    if args.abstract_product_counterexample:
        run_abstract_counterexample(); did = True
    if args.small_trees is not None:
        run_small_trees(args.small_trees); did = True
    if not did:
        ap.print_help()


if __name__ == "__main__":
    main()
