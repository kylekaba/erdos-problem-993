# Erdos Problem 993 (Open)

This repository collects exact-arithmetic computations and proof-route notes from two working sessions on [Erdos Problem #993](https://www.erdosproblems.com/993): whether the independent-set sequence of every tree or forest is unimodal.

The project is still in progress. As of July 2026, the official Erdős Problems page lists #993 as open and "falsifiable", meaning that one finite tree or forest with a non-unimodal independent-set sequence would disprove it. The scripts here instead found and verified a useful obstruction to one tempting proof strategy, plus a more promising induction framework to investigate next.

## Main session findings

Let

```text
I(G; x) = sum_k i_k(G) x^k
mu_k(G) = (k + 1) i_{k+1}(G) / i_k(G).
```

A natural attempted smoothness lemma,

```text
mu_{k+1} <= mu_k + 1,
```

is false for explicit trees.

The clean witness family is the bush `B(c,m)`: one hub, `m` branches adjacent to the hub, and each branch carrying `c` pendant paths of length 2. Its independence polynomial is

```text
I(B(c,m); x) = ((1 + 2x)^c + x(1 + x)^c)^m + x(1 + 2x)^(mc).
```

The smallest highlighted homogeneous witness is `B(8,14)`, with `n = 239`, independence number `alpha = 126`, and at `k = 113`,

```text
mu_114 - mu_113 = 1.082785... > 1.
```

Larger witnesses in `erdos993_mu_increment_counterexamples.json`, including `B(8,15)`, `B(10,36)`, and `B(12,80)`, push the increment higher. These witnesses remain unimodal and have no hinge violations. The obstruction is therefore not a counterexample to Erdős 993; it only rules out proof approaches that rely on bounded local smoothness of `mu_k`.

The apparent mechanism is a hub phase transition: independent sets containing the hub have a hard support ceiling, and when that phase disappears between adjacent levels, `mu_k` can jump sharply while the coefficient sequence stays unimodal.

## More promising proof route

The better target appears to be the ratio

```text
q_k(G) = i_{k+1}(G) / i_k(G),
```

with the goal of proving a no-upward-crossing principle: once `q_k <= 1`, later ratios should not cross back above `1`.

The proposed induction is based on the root decomposition

```text
I(F; x) = I(F - r; x) + x I(F - N[r]; x).
```

The compact root-deletion target is:

```text
I(F-r)_k < I(F-r)_{k-1}
    => I(F)_{k+1} <= I(F)_k.                      (RD)
```

By contraposition this gives single-root interlacing:

```text
I(F)_{k+1} > I(F)_k
    => I(F-r)_k >= I(F-r)_{k-1}.                  (SRI)
```

For a leaf `ell` of a tree `T`, with neighbor `u`, put

```text
A = T - ell
H = T - {ell, u}
J = T - N[u].
```

Then

```text
I(A; x) = I(H; x) + x I(J; x)
I(T; x) = I(A; x) + x I(H; x).
```

If a strict descent in `T` forces both `A` and `H` to be non-rising at the same level, induction would prove unimodality.

## Two missing lemmas

The current proof route reduces to two forest-specific inequalities.

1. **Admissible star-product IRD closure (ASPC).**

   For sequences `A = (a_j)`, `B = (b_j)`, define `C = A + xB`, so `c_j = a_j + b_{j-1}`. Interval right-dominance (IRD) is:

   ```text
   a_t < a_{s-1}
       => a_{t+1} + b_t <= a_s + b_{s-1}
   ```

   for every `1 <= s <= t`. The local case `s = t = k` is RD.

   Rooted-tree pairs are generated recursively from `(1+x, 1)` by

   ```text
   (P_i, Q_i) -> (prod_i P_i + x prod_i Q_i, prod_i P_i).
   ```

   The needed theorem is that if `(P_i, Q_i)` are admissible rooted-tree pairs, then `(prod_i P_i, prod_i Q_i)` satisfies IRD. Plain arbitrary product closure is false; the included `abstract_product_closure_counterexample.py` verifies an explicit sequence-level counterexample.

2. **Two-step no-hidden-rise (2NHR).**

   For every rooted forest pair `(H,J) = (F-r, F-N[r])`, prove a statement of the form:

   ```text
   h_{k+1} > h_k >= h_{k-1}
       => h_{k+1} + j_k > h_{k-1} + j_{k-1}.
   ```

   This would prevent the second smaller forest in the leaf transfer from still rising while the two-step extension has already begun to fall.

If ASPC and 2NHR are proved, the induction should be short: derive leaf transfer under strict descent, then induct on the number of vertices.

## Repository contents

```text
erdos993_mu_increment_counterexamples.json
    Exact coefficient data for homogeneous bush witnesses.

verify_993_counterexample.py
    Original standalone dependency-free verifier for B(8,14).

erdos993_handoff/HANDOFF.md
    Detailed research memo from the sessions.

erdos993_handoff/README.md
    Quick package-level run instructions.

erdos993_handoff/erdos993_common.py
    Shared exact polynomial, forest, bush, unimodality, and mu utilities.

erdos993_handoff/erdos993_tools.py
    Expanded utilities for exact forest computations and lemma checks.

erdos993_handoff/verify_bush_counterexamples.py
    Exact verifier for the uploaded JSON witness records.

erdos993_handoff/run_handoff_checks.py
    Command-line runner for the main checks.

erdos993_handoff/abstract_product_closure_counterexample.py
    Verifies that naive abstract RD product closure is false.

erdos993_handoff/lemma_lab.py
    Independent small-tree testbed for RD, SRI, IRD, leaf transfer, and 2NHR.
```

## Running checks

The bush and JSON checks use only Python's standard library:

```bash
python3 verify_993_counterexample.py
cd erdos993_handoff
python3 verify_bush_counterexamples.py --json ../erdos993_mu_increment_counterexamples.json
python3 abstract_product_closure_counterexample.py
```

The exhaustive non-isomorphic tree checks require NetworkX:

```bash
python3 -m pip install -r requirements.txt
cd erdos993_handoff
python3 run_handoff_checks.py --all
python3 lemma_lab.py --max-n 12
```

`run_handoff_checks.py --all` verifies the bush witnesses, checks the uploaded JSON data, displays the abstract product-closure counterexample, and tests the proposed lemmas on all non-isomorphic trees through 12 vertices.

## Future directions

The next expert should probably attack ASPC first. It is the bottleneck: ordinary convolution closure fails, but every rooted-tree pair has extra recursive structure that may carry a stronger sign-regularity or interval-majorization certificate.

Promising routes:

- A Cauchy-Binet or sign-regularity proof in which each rooted tree carries a two-row certificate stronger than IRD and preserved by admissible star products.
- A discrete matching or Morse-theoretic injection on the independence complex of a forest, root-compatible enough to imply IRD and 2NHR.
- A hard-core model argument, but one that proves lagged interval dominance rather than pointwise monotone likelihood ratio.

Less promising routes:

- Do not try to prove log-concavity for all trees; known counterexamples to tree independence-sequence log-concavity exist.
- Do not try to rescue the `mu_{k+1} <= mu_k + C` idea as the main proof engine. The bush family suggests the defect comes from global phase collapse, not a small local error term.

## References

- [Erdos Problem #993](https://www.erdosproblems.com/993)
- Abdul Basit and David Galvin, [On the independent set sequence of a tree](https://arxiv.org/abs/2006.12562)
- David Galvin and Justin Hilyard, [The independent set sequence of some families of trees](https://arxiv.org/abs/1701.02204)
