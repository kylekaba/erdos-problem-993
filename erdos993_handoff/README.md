# Erdős 993 handoff verification package

This package contains exact-arithmetic scripts supporting the handoff memo.

## Quick start

```bash
cd erdos993_handoff
python run_handoff_checks.py --all
```

The `--all` suite verifies the bush witnesses, checks the uploaded JSON data if it is one directory above the package, displays the abstract product-closure counterexample, and tests the proposed lemmas on all non-isomorphic trees through 12 vertices.

## Individual checks

```bash
python run_handoff_checks.py --bush
python run_handoff_checks.py --json ../erdos993_mu_increment_counterexamples.json
python run_handoff_checks.py --abstract-product-counterexample
python run_handoff_checks.py --small-trees 14
```

The small-tree enumeration uses NetworkX.  All arithmetic in the polynomial computations is exact Python integer/rational arithmetic.

## Main concepts tested

- Bush counterexamples to `mu_{k+1} <= mu_k + 1`.
- Unimodality and hinge checks for those bushes.
- The compact root-deletion dominance target:
  `I(F-r)_k < I(F-r)_{k-1} => I(F)_{k+1} <= I(F)_k`.
- The single-root interlacing target:
  `I(F)_{k+1} > I(F)_k => I(F-r)_k >= I(F-r)_{k-1}`.
- Root-pair interval right-dominance.
- Strict leaf no-hidden-rise tests.
- A sequence-level counterexample showing RD is not closed under arbitrary products.
