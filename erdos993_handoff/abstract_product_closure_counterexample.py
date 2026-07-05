#!/usr/bin/env python3
"""Exact counterexample to naive abstract product closure of RD.

RD(A,B): a_k < a_{k-1} => a_{k+1}+b_k <= a_k+b_{k-1}.
Even if two pairs satisfy RD, their product pair need not.
This demonstrates why a proof for Erdős 993 must use forest-admissible structure,
not arbitrary convolution closure.
"""
from __future__ import annotations

from erdos993_common import coeff, is_log_concave, poly_mul


def rd_holds(A, B):
    for k in range(1, max(len(A), len(B)) + 2):
        if coeff(A, k) < coeff(A, k - 1):
            if coeff(A, k + 1) + coeff(B, k) > coeff(A, k) + coeff(B, k - 1):
                return False, k
    return True, None


def main() -> None:
    A = [1, 7, 8, 8, 8, 7]
    B = [1, 2, 4, 8, 8, 7]
    C = [1, 6, 8]
    D = [1, 4, 7]
    AC = poly_mul(A, C)
    BD = poly_mul(B, D)

    for name, seq in [("A", A), ("B", B), ("C", C), ("D", D), ("AC", AC), ("BD", BD)]:
        print(f"{name} = {seq}; log_concave={is_log_concave(seq)}")

    for name, X, Y in [("(A,B)", A, B), ("(C,D)", C, D), ("(AC,BD)", AC, BD)]:
        ok, k = rd_holds(X, Y)
        print(f"RD{name}: {ok}" + (f"; first failure k={k}" if not ok else ""))

    k = 5
    print("\nProduct failure at k=5:")
    print(f"AC_5={coeff(AC,5)} < AC_4={coeff(AC,4)}")
    print(f"AC_6+BD_5={coeff(AC,6)}+{coeff(BD,5)}={coeff(AC,6)+coeff(BD,5)}")
    print(f"AC_5+BD_4={coeff(AC,5)}+{coeff(BD,4)}={coeff(AC,5)+coeff(BD,4)}")


if __name__ == "__main__":
    main()
