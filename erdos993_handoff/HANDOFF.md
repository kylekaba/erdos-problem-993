# Erdős Problem 993 handoff memo

## Problem

Erdős problem 993 asks whether the independent-set sequence of every tree or forest is unimodal.  If

\[
I(G;x)=\sum_k i_k(G)x^k,
\]

then the question is whether every tree/forest has

\[
i_0\le i_1\le\cdots\le i_m\ge i_{m+1}\ge\cdots.
\]

As of this handoff, the problem is still listed as open on erdosproblems.com.

## Main negative result from this session

A tempting strengthening fails badly.  Define

\[
\mu_k(G)=\frac{(k+1)i_{k+1}(G)}{i_k(G)}.
\]

The hoped-for lemma

\[
\mu_{k+1}\le \mu_k+1
\]

is false for explicit trees.

The clean witness is the bush \(B(c,m)\): one hub vertex, \(m\) adjacent branch vertices, and each branch vertex supports \(c\) pendant paths of length 2.  Its independence polynomial is

\[
I(B(c,m);x)=\left((1+2x)^c+x(1+x)^c\right)^m+x(1+2x)^{mc}.
\]

For \(B(8,14)\), \(n=239\), \(\alpha=126\), and at \(k=113\),

\[
\mu_{114}-\mu_{113}=1.082785\ldots>1.
\]

The same witness remains unimodal and has no hinge violations.  Larger bushes give larger increments, e.g. \(B(12,80)\) has increment about \(2.404711\).  This kills any proof route that depends on bounded smoothness of \(\mu_k\).

## Correct lesson from the bush obstruction

The failure is a phase-coexistence effect.  The hub-in phase has a sharp support ceiling.  At one level, many independent sets include the hub and have zero availability; at the next level, the hub-in phase disappears.  This can make \(\mu_k\) jump even though unimodality is intact.

Therefore the target should not be a smooth increment bound for \(\mu\).  The target should be a no-upward-crossing principle for the ratios

\[
q_k(G)=\frac{i_{k+1}(G)}{i_k(G)}.
\]

Unimodality is equivalent to saying that after \(q_k\) drops below 1, it never rises above 1 again; for a proof by induction it is enough to show that after a strict descent in the coefficient sequence, all later steps remain nonincreasing.

## Main positive proof architecture

Let \(r\) be a vertex of a forest \(F\).  Write

\[
I(F;x)=I(F-r;x)+xI(F-N[r];x).
\]

The compact target is the root-deletion dominance statement:

\[
I(F-r)_k<I(F-r)_{k-1}\quad\Longrightarrow\quad I(F)_{k+1}\le I(F)_k.
\tag{RD}
\]

If RD is true, it immediately implies the single-root interlacing lemma:

\[
I(F)_{k+1}>I(F)_k\quad\Longrightarrow\quad I(F-r)_k\ge I(F-r)_{k-1}.
\tag{SRI}
\]

Indeed, the contrapositive of RD says that if \(I(F)_{k+1}>I(F)_k\), then \(I(F-r)_k\not<I(F-r)_{k-1}\).

## Leaf-transfer route to the original problem

Let \(\ell\) be a leaf of \(T\), and \(u\) its neighbor.  Put

\[
A=T-\ell,\qquad H=T-\{\ell,u\},\qquad J=T-N[u].
\]

Then

\[
I(A;x)=I(H;x)+xI(J;x),
\]

and

\[
I(T;x)=I(A;x)+xI(H;x).
\]

A strict descent in \(T\) at level \(k\) should force strict/no-rise in both smaller forests:

\[
I(T)_{k+1}<I(T)_k\quad\Longrightarrow\quad I(A)_{k+1}\le I(A)_k,
\]

and

\[
I(T)_{k+1}<I(T)_k\quad\Longrightarrow\quad I(H)_{k+1}\le I(H)_k.
\]

The first implication follows from RD/SRI.  The second needs a companion two-step no-hidden-rise lemma.

If both implications hold, induction proves unimodality: once \(T\) strictly descends, both smaller forests are at or beyond their modes by induction, so \(T\) cannot rise again.

## Interval right-dominance

For sequences \(A=(a_k)\), \(B=(b_k)\), define \(C=A+xB\), i.e.

\[
c_k=a_k+b_{k-1}.
\]

The proposed strengthening is interval right-dominance (IRD): for all \(1\le s\le t\),

\[
a_t<a_{s-1}\quad\Longrightarrow\quad a_{t+1}+b_t\le a_s+b_{s-1}.
\tag{IRD}
\]

The local case \(s=t=k\) is RD.

The base pair \((1+x,1)\) satisfies IRD by direct inspection.

## The sharp missing theorem

Plain product closure for RD/IRD is false.  An explicit abstract counterexample is built into the script:

\[
A=(1,7,8,8,8,7),\quad B=(1,2,4,8,8,7),
\]

\[
C=(1,6,8),\quad D=(1,4,7).
\]

Both pairs satisfy RD, but \((AC,BD)\) fails RD.

Therefore the missing theorem must use forest-admissibility, not arbitrary convolution.

Rooted tree pairs are generated recursively by the operation

\[
(P_i,Q_i)_{i=1}^m\mapsto\left(\prod_iP_i+x\prod_iQ_i,\ \prod_iP_i\right),
\]

starting from \((1+x,1)\).  The needed forest-specific theorem is:

> If \((P_i,Q_i)\) are admissible rooted-tree pairs, then \((\prod_i P_i,\prod_i Q_i)\) satisfies IRD.

Call this admissible star-product closure (ASPC).  ASPC implies RD for every rooted forest, hence SRI.

## Path/spider evidence

For endpoint-rooted paths, the coefficients are

\[
I(P_n;x)=\sum_k \binom{n-k+1}{k}x^k.
\]

A direct ratio calculation proves RD for endpoint-rooted path extension.  The kernel

\[
K(n,k)=\binom{n-k+1}{k}
\]

is the total-positivity object that should underwrite the spider case by variation-diminishing convolution.  Thus the proposed route works cleanly for paths and likely for spiders.  General trees require a branched/recursive total-positivity certificate.

## Cleanest path to a full proof

1. Prove ASPC: forest-admissible star products preserve IRD.
2. Prove the two-step no-hidden-rise lemma:
   \[
   h_{k+1}>h_k\ge h_{k-1}\quad\Longrightarrow\quad h_{k+1}+j_k>h_{k-1}+j_{k-1}
   \]
   for every rooted forest pair \((H,J)=(F-r,F-N[r])\), or find an equivalent IRD consequence that gives the same leaf-transfer conclusion.
3. Deduce leaf-transfer under strict descent.
4. Run the standard induction on vertex count.

This would prove Erdős 993 for forests directly, without relying on log-concavity or smoothness of \(\mu_k\).

## Files in this package

- `erdos993_tools.py`: exact integer polynomial utilities, bush formulas, RD/IRD tests, root/leaf lemma tests.
- `run_handoff_checks.py`: command-line runner for the main verification tasks.
- `../erdos993_mu_increment_counterexamples.json`: uploaded exact witness data.
- `../verify_993_counterexample.py`: original standalone verifier for the clean bush witness.

Suggested checks:

```bash
python run_handoff_checks.py --bush
python run_handoff_checks.py --json ../erdos993_mu_increment_counterexamples.json
python run_handoff_checks.py --abstract-product-counterexample
python run_handoff_checks.py --small-trees 12
```
