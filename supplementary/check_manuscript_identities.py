#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Selected algebra checks for the merged Johnson core manuscript.

The script uses exact SymPy expressions for formal identities and exact
integer arithmetic for finite-range regression checks.  It does not
mechanize the combinatorial or representation-theoretic reductions.
"""

from __future__ import annotations

from math import comb

import sympy as sp


def require_zero(name: str, expr: sp.Expr) -> None:
    reduced = sp.factor(sp.cancel(expr))
    if reduced != 0:
        raise AssertionError(f"{name} failed: {reduced}")
    print(f"PASS  {name}")


def require_equal(name: str, left: sp.Expr, right: sp.Expr) -> None:
    require_zero(name, left - right)


def require(condition: bool, name: str, detail: object | None = None) -> None:
    if not condition:
        suffix = "" if detail is None else f": {detail}"
        raise AssertionError(f"{name} failed{suffix}")


def eberlein(k: int, j: int, i: int) -> int:
    total = 0
    for r in range(i + 1):
        # Guard invalid arguments because math.comb requires 0 <= r <= n.
        if not (0 <= i - r <= k - r):
            continue
        if not (0 <= r <= k - j):
            continue
        total += (
            (-1) ** (i - r)
            * comb(k - r, i - r)
            * comb(k - j, r)
            * comb(k + 1 + r - j, r)
        )
    return total


def check_endpoint_eigenvalues() -> None:
    for k in range(4, 41):
        K = k + 1
        for j in range(k + 1):
            t = K - j
            eps = (-1) ** j
            expected = {
                1: t * t - K,
                k - 1: eps * t * (t * t - 2 * K + 1) // 2,
                k: eps * t,
            }
            for i, value in expected.items():
                actual = eberlein(k, j, i)
                if actual != value:
                    raise AssertionError(
                        f"endpoint formula failed at k={k}, j={j}, i={i}: "
                        f"{actual} != {value}"
                    )
    print("PASS  endpoint Eberlein formulas (k=4..40, all j)")


K, t, s = sp.symbols("K t s")
X, Y = sp.symbols("X Y")
Z = 1 - X - Y


def general_e(eps: int, u: sp.Expr, fibre_size: sp.Expr = s) -> sp.Expr:
    a = (u**2 - K) / (K * (K - 1))
    b = eps * u * (u**2 - 2 * K + 1) / (K * (K - 1) ** 2)
    c = eps * u / K
    return sp.factor(1 + (fibre_size - 1) * (a * X + b * Y + c * Z))


def specialized_e(eps: int, u: sp.Expr) -> sp.Expr:
    return sp.factor(
        1
        + 2 * (u**2 - K) / (K - 1) * X
        + 2 * eps * u * (u**2 - 2 * K + 1) / (K - 1) ** 2 * Y
        + 2 * eps * u * Z
    )


def check_full_case_specialization() -> None:
    fibre_size = 2 * K + 1
    require_equal(
        "[eq:full-bracket] full-case specialization, negative parity",
        specialized_e(-1, t),
        general_e(-1, t, fibre_size),
    )
    require_equal(
        "[eq:full-bracket] full-case specialization, positive parity",
        specialized_e(+1, t),
        general_e(+1, t, fibre_size),
    )


def check_proper_divisor_identities() -> None:
    a_plus = (t**2 - 1) / (K * (K - 1))
    b_plus = (t**3 - (2 * K - 1) * t + (K - 1) ** 2) / (
        K * (K - 1) ** 2
    )
    c_plus = (t + 1) / K
    even_rhs = (K - s + 1) / K + (s - 1) * (
        a_plus * X + b_plus * Y + c_plus * Z
    )
    require_equal(
        "[eq:T1-even-identity] proper-divisor even-module identity",
        general_e(+1, t),
        even_rhs,
    )

    A_t = (t - 1) * (t + K) / (K * (K - 1))
    C_t = t * (K**2 - t**2) / (K * (K - 1) ** 2)
    D = K - (s - 1) * t
    odd_rhs = D / K + (s - 1) * (A_t * X + C_t * Y)
    require_equal(
        "[eq:T1-odd-basic] proper-divisor odd-module identity",
        general_e(-1, t),
        odd_rhs,
    )

    # The D=0 identity is evaluated on s-1=K/t and j=K-t.
    s_d0 = 1 + K / t
    j_minus_one = K - t - 1
    lambda_t = (t - 1) * (t + K) / ((K - 2) * (2 * K - 1))
    coefficient = (
        (K + t)
        * (K - t - 1)
        * (t * (K - 2) + 1)
        / (t * (K - 2) * (K - 1) ** 2)
    )
    d0_rhs = (
        lambda_t * general_e(-1, K - 1, s_d0)
        + lambda_t * j_minus_one / t
        + coefficient * Y
    )
    require_equal(
        "[eq:T1-D0] D_j=0 auxiliary affine identity",
        general_e(-1, t, s_d0),
        d0_rhs,
    )


def check_odd_K_interpolation() -> None:
    B1 = specialized_e(-1, K - 1)
    Bkm1 = specialized_e(-1, 2)
    Bk = specialized_e(+1, 1)

    am = (t - 2) * (t + 1) * (K * (K - t - 1) + 4 * t - 3) / (
        K * (K - 3) * (4 * K - 7)
    )
    bm = (t + 1) * (K - t - 1) * (t * (2 * K - 5) + 3) / (
        3 * (K - 3) * (4 * K - 7)
    )
    gm = (t - 2) * (K - t - 1) * (4 * t * (K - 1) + 3) / (
        3 * K * (4 * K - 7)
    )
    require_equal(
        "[eq:odd-minus-interp] endpoint interpolation, odd K, negative parity",
        specialized_e(-1, t),
        am * B1 + bm * Bkm1 + gm * Bk,
    )

    ap = (t - 1) * (t + 2) * (K * (K - 1) + t * (K - 4) - 3) / (
        K * (K - 3) * (4 * K - 7)
    )
    bp = (t - 1) * (K + t - 1) * (t * (2 * K - 5) - 3) / (
        3 * (K - 3) * (4 * K - 7)
    )
    gp = (t + 2) * (K + t - 1) * (4 * t * (K - 1) - 3) / (
        3 * K * (4 * K - 7)
    )
    require_equal(
        "[eq:odd-plus-interp] endpoint interpolation, odd K, positive parity",
        specialized_e(+1, t),
        ap * B1 + bp * Bkm1 + gp * Bk,
    )

def check_even_K_interpolation() -> None:
    B1 = specialized_e(-1, K - 1)
    Bkm1 = specialized_e(+1, 2)
    Bk = specialized_e(-1, 1)

    am = (t - 1) * (t + 2) * (K * (K - t + 1) + 2 * t - 1) / (
        (K - 2) * (K + 1) * (4 * K - 3)
    )
    bm = (t - 1) * (K - t - 1) * (t * (2 * K - 1) + 1) / (
        3 * (K + 1) * (4 * K - 3)
    )
    gm = (t + 2) * (K - t - 1) * (4 * t * (K - 2) + 5) / (
        3 * (K - 2) * (4 * K - 3)
    )
    require_equal(
        "[eq:even-minus-interp] endpoint interpolation, even K, negative parity",
        specialized_e(-1, t),
        am * B1 + bm * Bkm1 + gm * Bk,
    )

    ap = (t - 2) * (t + 1) * (K * (K + t + 1) - 2 * t - 1) / (
        (K - 2) * (K + 1) * (4 * K - 3)
    )
    bp = (t + 1) * (K + t - 1) * (t * (2 * K - 1) - 1) / (
        3 * (K + 1) * (4 * K - 3)
    )
    gp = (t - 2) * (K + t - 1) * (4 * t * (K - 2) - 5) / (
        3 * (K - 2) * (4 * K - 3)
    )
    require_equal(
        "[eq:even-plus-interp] endpoint interpolation, even K, positive parity",
        specialized_e(+1, t),
        ap * B1 + bp * Bkm1 + gp * Bk,
    )

def check_core_identities() -> None:
    q = sp.symbols("q")
    fibre_size = 2 * K + 1
    u = sp.symbols("u")
    total_nonconstant = q * (fibre_size - 1) / fibre_size
    v = total_nonconstant - u
    h0 = q / fibre_size

    # Odd K: support j=k-1 (epsilon=-1,t=2) and j=k (epsilon=+1,t=1).
    R1_odd = K * (K - 1) * h0 + (4 - K) * u + (1 - K) * v
    Rk_odd = K * h0 - 2 * u + v
    require_equal(
        "[eq:odd-core-identity] odd-K endpoint core identity",
        R1_odd + Rk_odd,
        q / fibre_size * K * (4 - K),
    )

    # Even K: support j=k-1 (epsilon=+1,t=2) and j=k (epsilon=-1,t=1).
    R1_even = K * (K - 1) * h0 + (4 - K) * u + (1 - K) * v
    require_equal(
        "[eq:even-core-identity] even-K endpoint core identity",
        R1_even + 3 * v,
        q / fibre_size * K * (7 - K),
    )


def interpolation_numerators(kval: int, tv: int, eps: int) -> tuple[int, int, int]:
    if kval % 2 and eps == -1:
        return (
            (tv - 2) * (tv + 1) * (kval * (kval - tv - 1) + 4 * tv - 3),
            (tv + 1) * (kval - tv - 1) * (tv * (2 * kval - 5) + 3),
            (tv - 2) * (kval - tv - 1) * (4 * tv * (kval - 1) + 3),
        )
    if kval % 2 and eps == 1:
        return (
            (tv - 1)
            * (tv + 2)
            * (kval * (kval - 1) + tv * (kval - 4) - 3),
            (tv - 1) * (kval + tv - 1) * (tv * (2 * kval - 5) - 3),
            (tv + 2) * (kval + tv - 1) * (4 * tv * (kval - 1) - 3),
        )
    if kval % 2 == 0 and eps == -1:
        return (
            (tv - 1)
            * (tv + 2)
            * (kval * (kval - tv + 1) + 2 * tv - 1),
            (tv - 1) * (kval - tv - 1) * (tv * (2 * kval - 1) + 1),
            (tv + 2) * (kval - tv - 1) * (4 * tv * (kval - 2) + 5),
        )
    if kval % 2 == 0 and eps == 1:
        return (
            (tv - 2)
            * (tv + 1)
            * (kval * (kval + tv + 1) - 2 * tv - 1),
            (tv + 1) * (kval + tv - 1) * (tv * (2 * kval - 1) - 1),
            (tv - 2) * (kval + tv - 1) * (4 * tv * (kval - 2) - 5),
        )
    raise AssertionError(f"unexpected parity data: K={kval}, t={tv}, eps={eps}")


def interpolation_denominators(kval: int) -> tuple[int, int, int]:
    if kval % 2:
        return (
            kval * (kval - 3) * (4 * kval - 7),
            3 * (kval - 3) * (4 * kval - 7),
            3 * kval * (4 * kval - 7),
        )
    return (
        (kval - 2) * (kval + 1) * (4 * kval - 3),
        3 * (kval + 1) * (4 * kval - 3),
        3 * (kval - 2) * (4 * kval - 3),
    )


def check_coefficient_grids() -> None:
    """Finite regression check of parity grids and coefficient signs.

    Universal signs are proved by the factorizations in Appendix A.
    """

    for kval in range(5, 202):
        derived = {((-1) ** j, kval - j) for j in range(1, kval)}
        if kval % 2:
            expected = {
                *{(-1, tv) for tv in range(2, kval, 2)},
                *{(+1, tv) for tv in range(1, kval - 1, 2)},
            }
            endpoints = {(-1, kval - 1), (-1, 2), (+1, 1)}
        else:
            expected = {
                *{(-1, tv) for tv in range(1, kval, 2)},
                *{(+1, tv) for tv in range(2, kval, 2)},
            }
            endpoints = {(-1, kval - 1), (+1, 2), (-1, 1)}

        require(
            derived == expected,
            "parity-grid derivation",
            {"K": kval, "derived": derived, "expected": expected},
        )
        denominators = interpolation_denominators(kval)
        require(
            all(value > 0 for value in denominators),
            "interpolation denominator positivity",
            {"K": kval, "denominators": denominators},
        )
        for eps, tv in derived:
            values = interpolation_numerators(kval, tv, eps)
            require(
                all(value >= 0 for value in values),
                "coefficient nonnegativity regression",
                {"K": kval, "t": tv, "eps": eps, "numerators": values},
            )
            if (eps, tv) not in endpoints:
                require(
                    values[0] > 0,
                    "strict B_1 coefficient regression",
                    {"K": kval, "t": tv, "eps": eps, "alpha": values[0]},
                )
    print("PASS  parity/sign regression checks (K=5..201)")


def main() -> None:
    check_endpoint_eigenvalues()
    check_full_case_specialization()
    check_proper_divisor_identities()
    check_odd_K_interpolation()
    check_even_K_interpolation()
    check_core_identities()
    check_coefficient_grids()
    print("All checks in check_manuscript_identities.py completed successfully.")


if __name__ == "__main__":
    main()
