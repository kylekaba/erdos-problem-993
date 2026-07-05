#!/usr/bin/env python3
"""Standalone, dependency-free verification that the inequality
        mu_{k+1} <= mu_k + 1,   mu_k = (k+1) i_{k+1} / i_k
(equivalently  Var|F-N[S]| <= 2 E|F-N[S]| + 2 E|E(F-N[S])|  over uniform S in I_k(F))
FAILS for an explicit tree: the "bush" B(c,m) with c=8, m=14 (n=239 vertices).

B(c,m): hub v0 adjacent to v_1..v_m; each v_i additionally adjacent to c pendant
paths of length 2 (v_i - a - b).  n = m(2c+1) + 1.

Two independent computations of the i-sequence:
  (1) exact integer DP over the tree (works for any tree), and
  (2) the closed form  I(B) = [(2x+1)^c + x(1+x)^c]^m + x(2x+1)^{mc},
      from deleting/contracting the hub:  I(G) = I(G - v0) + x I(G - N[v0]).
All arithmetic is exact (Python ints / Fractions).
"""
from fractions import Fraction

C, M = 8, 14

# ---------- build the tree ----------
edges = []
nid = 1  # 0 is the hub
for _ in range(M):
    vi = nid; nid += 1
    edges.append((0, vi))
    for _ in range(C):
        a, b = nid, nid + 1; nid += 2
        edges += [(vi, a), (a, b)]
N = nid
adj = [[] for _ in range(N)]
for u, v in edges:
    adj[u].append(v); adj[v].append(u)
assert N == M * (2 * C + 1) + 1

# ---------- (1) exact DP ----------
def indep_seq(n, adj, root=0):
    parent = {root: -1}; order = [root]; seen = [False] * n; seen[root] = True
    i = 0
    while i < len(order):
        v = order[i]; i += 1
        for u in adj[v]:
            if not seen[u]:
                seen[u] = True; parent[u] = v; order.append(u)
    def pm(a, b):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    r[i + j] += x * y
        return r
    def pa(a, b):
        if len(a) < len(b): a, b = b, a
        return [a[i] + (b[i] if i < len(b) else 0) for i in range(len(a))]
    p0 = {v: [1] for v in order}; p1 = {v: [0, 1] for v in order}
    for v in reversed(order):
        for u in adj[v]:
            if u != parent[v]:
                p1[v] = pm(p1[v], p0[u])
                p0[v] = pm(p0[v], pa(p0[u], p1[u]))
    return pa(p0[root], p1[root])

s_dp = indep_seq(N, adj)

# ---------- (2) closed form ----------
def pm(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r
def pa(a, b):
    if len(a) < len(b): a, b = b, a
    return [a[i] + (b[i] if i < len(b) else 0) for i in range(len(a))]
def pp(p, e):
    r = [1]
    while e:
        if e & 1: r = pm(r, p)
        p = pm(p, p); e >>= 1
    return r
arm = pa(pp([1, 2], C), pm([0, 1], pp([1, 1], C)))       # (2x+1)^c + x(1+x)^c
s_cf = pa(pp(arm, M), pm([0, 1], pp([1, 2], M * C)))     # (arm)^m + x(2x+1)^{mc}

assert s_dp == s_cf, "DP and closed form disagree!"
s = s_dp
alpha = len(s) - 1

# ---------- the violation ----------
mu = [Fraction((k + 1) * s[k + 1], s[k]) for k in range(alpha)]
k = max(range(alpha - 1), key=lambda k: mu[k + 1] - mu[k])
inc = mu[k + 1] - mu[k]
print(f"tree: bush(c={C}, m={M}), n = {N}, alpha = {alpha}")
print(f"k = {k}:   mu_k     = {float(mu[k]):.6f}")
print(f"          mu_(k+1) = {float(mu[k+1]):.6f}")
print(f"          increment = {float(inc):.6f}   exact fraction > 1: {inc > 1}")
print(f"variance-form defect  VarX - 2EX - 2Ee = mu_k*(inc - 1) = {float(mu[k]*(inc-1)):.6f}  (> 0)")

# unimodality of the same tree (sanity: unimodality itself is NOT violated)
drop = next(j for j in range(alpha) if s[j + 1] < s[j])
assert all(s[j + 1] <= s[j] for j in range(drop, alpha)), "non-unimodal?!"
print("i-sequence of this tree IS unimodal; the hinge implication "
      "(mu_k <= k+1  =>  mu_(k+1) <= k+2) also holds at every k.")
hv = [j for j in range(alpha - 1) if mu[j] <= j + 1 and mu[j + 1] > j + 2]
assert not hv
assert inc > 1
print("VERIFIED: mu_{k+1} <= mu_k + 1 is FALSE for this forest.")
