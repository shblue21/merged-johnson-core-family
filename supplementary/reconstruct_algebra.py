#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Separate algebra reconstruction checks for the core manuscript.

The program reconstructs the interpolation coefficients for the case
s=2k+3 by solving the defining affine 3-by-3 systems.  It then compares
those separately derived answers with the factored expressions in the
paper.  It also performs finite-range endpoint, arithmetic, and sign
regression checks.

Its scope is algebraic: it does not mechanize Schur orthogonality,
the minimum-core reduction, or the extremal set-system proof.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, gcd

import sympy as sp


def assert_zero(label: str, expression: sp.Expr) -> None:
    value = sp.factor(sp.cancel(expression))
    if value != 0:
        raise AssertionError(f"{label}: {value}")
    print(f"PASS  {label}")


def require(condition: bool, label: str, detail: object | None = None) -> None:
    if not condition:
        suffix = "" if detail is None else f": {detail}"
        raise AssertionError(f"{label} failed{suffix}")


def eberlein(k: int, j: int, i: int) -> int:
    """First-eigenmatrix entry from the defining finite sum."""
    total = 0
    for r in range(i + 1):
        if r > k - j:
            continue
        total += (
            (-1) ** (i - r)
            * comb(k - r, i - r)
            * comb(k - j, r)
            * comb(k + 1 + r - j, r)
        )
    return total


def endpoint_reconstruction() -> None:
    for k in range(4, 81):
        K0 = k + 1
        for j in range(k + 1):
            t0 = K0 - j
            eps = (-1) ** j
            reconstructed = (
                eberlein(k, j, 1),
                eberlein(k, j, k - 1),
                eberlein(k, j, k),
            )
            closed = (
                t0 * t0 - K0,
                eps * t0 * (t0 * t0 - 2 * K0 + 1) // 2,
                eps * t0,
            )
            if reconstructed != closed:
                raise AssertionError((k, j, reconstructed, closed))
        valencies = (
            eberlein(k, 0, 1),
            eberlein(k, 0, k - 1),
            eberlein(k, 0, k),
        )
        expected_valencies = (
            K0 * (K0 - 1),
            K0 * (K0 - 1) ** 2 // 2,
            K0,
        )
        require(
            valencies == expected_valencies,
            "endpoint valency reconstruction",
            {"k": k, "actual": valencies, "expected": expected_valencies},
        )
    print("PASS  Eberlein endpoint/valency reconstruction (k=4..80)")


K, t, X, Y = sp.symbols("K t X Y")
Z = 1 - X - Y


def specialized_bracket(eps: int, u: sp.Expr) -> sp.Expr:
    return sp.factor(
        1
        + 2 * (u**2 - K) / (K - 1) * X
        + 2 * eps * u * (u**2 - 2 * K + 1) / (K - 1) ** 2 * Y
        + 2 * eps * u * Z
    )


def affine_vector(form: sp.Expr) -> sp.Matrix:
    """Return coefficients of 1,X,Y after Z=1-X-Y."""
    poly = sp.Poly(sp.together(form), X, Y)
    return sp.Matrix(
        [
            poly.coeff_monomial(1),
            poly.coeff_monomial(X),
            poly.coeff_monomial(Y),
        ]
    )


def solve_coefficients(
    target: sp.Expr, left: sp.Expr, middle: sp.Expr, right: sp.Expr
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    basis = sp.Matrix.hstack(
        affine_vector(left), affine_vector(middle), affine_vector(right)
    )
    solution = basis.inv() * affine_vector(target)
    return tuple(sp.factor(sp.cancel(x)) for x in solution)


def declared_odd_minus() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    return (
        (t - 2) * (t + 1) * (K * (K - t - 1) + 4 * t - 3)
        / (K * (K - 3) * (4 * K - 7)),
        (t + 1) * (K - t - 1) * (t * (2 * K - 5) + 3)
        / (3 * (K - 3) * (4 * K - 7)),
        (t - 2) * (K - t - 1) * (4 * t * (K - 1) + 3)
        / (3 * K * (4 * K - 7)),
    )


def declared_odd_plus() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    return (
        (t - 1) * (t + 2) * (K * (K - 1) + t * (K - 4) - 3)
        / (K * (K - 3) * (4 * K - 7)),
        (t - 1) * (K + t - 1) * (t * (2 * K - 5) - 3)
        / (3 * (K - 3) * (4 * K - 7)),
        (t + 2) * (K + t - 1) * (4 * t * (K - 1) - 3)
        / (3 * K * (4 * K - 7)),
    )


def declared_even_minus() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    return (
        (t - 1) * (t + 2) * (K * (K - t + 1) + 2 * t - 1)
        / ((K - 2) * (K + 1) * (4 * K - 3)),
        (t - 1) * (K - t - 1) * (t * (2 * K - 1) + 1)
        / (3 * (K + 1) * (4 * K - 3)),
        (t + 2) * (K - t - 1) * (4 * t * (K - 2) + 5)
        / (3 * (K - 2) * (4 * K - 3)),
    )


def declared_even_plus() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    return (
        (t - 2) * (t + 1) * (K * (K + t + 1) - 2 * t - 1)
        / ((K - 2) * (K + 1) * (4 * K - 3)),
        (t + 1) * (K + t - 1) * (t * (2 * K - 1) - 1)
        / (3 * (K + 1) * (4 * K - 3)),
        (t - 2) * (K + t - 1) * (4 * t * (K - 2) - 5)
        / (3 * (K - 2) * (4 * K - 3)),
    )


def compare_triple(
    label: str,
    derived: tuple[sp.Expr, sp.Expr, sp.Expr],
    declared: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> None:
    require(
        len(derived) == len(declared) == 3,
        f"{label} coefficient-count",
        {"derived": len(derived), "declared": len(declared)},
    )
    for index, (left, right) in enumerate(zip(derived, declared), 1):
        assert_zero(f"{label} coefficient {index}", left - right)


def interpolation_reconstruction() -> None:
    odd_basis = (
        specialized_bracket(-1, K - 1),
        specialized_bracket(-1, 2),
        specialized_bracket(+1, 1),
    )
    compare_triple(
        "odd-K negative parity",
        solve_coefficients(specialized_bracket(-1, t), *odd_basis),
        declared_odd_minus(),
    )
    compare_triple(
        "odd-K positive parity",
        solve_coefficients(specialized_bracket(+1, t), *odd_basis),
        declared_odd_plus(),
    )

    even_basis = (
        specialized_bracket(-1, K - 1),
        specialized_bracket(+1, 2),
        specialized_bracket(-1, 1),
    )
    compare_triple(
        "even-K negative parity",
        solve_coefficients(specialized_bracket(-1, t), *even_basis),
        declared_even_minus(),
    )
    compare_triple(
        "even-K positive parity",
        solve_coefficients(specialized_bracket(+1, t), *even_basis),
        declared_even_plus(),
    )

def arithmetic_partition() -> None:
    full_cases = []
    proper_cases = []
    direct_cases = []
    for k in range(4, 501):
        K0 = k + 1
        m = 2 * k + 3
        N = comb(2 * k + 1, k)
        d = gcd(N, m)
        for s0 in range(2, d + 1):
            if d % s0:
                continue
            if s0 == m:
                full_cases.append((k, s0))
            else:
                require(s0 >= 3, "proper divisor lower bound", (k, s0))
                require(
                    s0 <= (2 * K0 + 1) // 3,
                    "proper divisor quotient bound",
                    (k, s0),
                )
                require(s0 <= K0, "proper divisor K bound", (k, s0))
                proper_cases.append((k, s0))
        for s0 in range(2, m + 1):
            if N % s0 == 0 and m % s0 == 0:
                direct_cases.append((k, s0))

    gcd_cases = sorted(proper_cases + full_cases)
    require(
        gcd_cases == sorted(direct_cases),
        "divisor enumeration cross-check",
        {"gcd_count": len(gcd_cases), "direct_count": len(direct_cases)},
    )
    require(
        (len(proper_cases), len(full_cases)) == (449, 13),
        "divisor enumeration coverage sentinel",
        {"proper": len(proper_cases), "full": len(full_cases)},
    )
    require((5, 13) not in full_cases, "K=6 full-case exclusion", full_cases)
    print(
        "PASS  divisor partition stress test "
        f"(k=4..500, proper={len(proper_cases)}, full={len(full_cases)})"
    )


def sign_domains() -> None:
    families = {
        "odd_minus": declared_odd_minus(),
        "odd_plus": declared_odd_plus(),
        "even_minus": declared_even_minus(),
        "even_plus": declared_even_plus(),
    }
    evaluators = {}
    for name, formulas in families.items():
        compiled = []
        for expression in formulas:
            numerator, denominator = sp.fraction(sp.cancel(expression))
            compiled.append(
                (
                    sp.lambdify((K, t), numerator, modules="math"),
                    sp.lambdify((K, t), denominator, modules="math"),
                )
            )
        evaluators[name] = compiled

    for K0 in range(5, 501):
        if K0 % 2:
            family_for_eps = {-1: "odd_minus", 1: "odd_plus"}
            endpoints = {(-1, K0 - 1), (-1, 2), (1, 1)}
        else:
            family_for_eps = {-1: "even_minus", 1: "even_plus"}
            endpoints = {(-1, K0 - 1), (1, 2), (-1, 1)}

        for j in range(1, K0):
            eps = (-1) ** j
            t0 = K0 - j
            compiled = evaluators[family_for_eps[eps]]
            values = []
            for numerator, denominator in compiled:
                denominator_value = denominator(K0, t0)
                require(
                    denominator_value > 0,
                    "interpolation denominator positivity",
                    {
                        "K": K0,
                        "j": j,
                        "t": t0,
                        "eps": eps,
                        "denominator": denominator_value,
                    },
                )
                values.append(Fraction(numerator(K0, t0), denominator_value))
            require(
                all(value >= 0 for value in values),
                "coefficient-sign regression",
                {"K": K0, "j": j, "t": t0, "eps": eps, "values": values},
            )
            if (eps, t0) not in endpoints:
                require(
                    values[0] > 0,
                    "strict B_1 coefficient regression",
                    {"K": K0, "j": j, "t": t0, "eps": eps, "alpha": values[0]},
                )
    print("PASS  coefficient-sign regression checks (K=5..500)")


def displayed_independent_family_construction() -> None:
    for k in range(4, 21):
        U = set(range(k + 1))
        T = set(range(k + 1, 2 * k))
        p = 2 * k
        family = []
        for missing in U:
            family.append(frozenset(U - {missing}))
        for x in U | {p}:
            family.append(frozenset(T | {x}))
        require(
            len(family) == 2 * k + 3,
            "displayed family size",
            {"k": k, "size": len(family)},
        )
        require(
            len(set(family)) == len(family),
            "displayed family distinctness",
            {"k": k, "size": len(family), "distinct": len(set(family))},
        )
        for i, left in enumerate(family):
            for right in family[i + 1 :]:
                intersection = len(left & right)
                require(
                    intersection in {0, 1, k - 1},
                    "displayed family intersection condition",
                    {"k": k, "intersection": intersection, "left": left, "right": right},
                )
    print("PASS  displayed independent-family construction (k=4..20)")


def core_endpoint_identities_from_eberlein() -> None:
    for K0 in range(5, 151):
        k = K0 - 1
        s0 = 2 * K0 + 1
        q = Fraction(1, 1)
        h0 = q / s0
        total = q * (s0 - 1) / s0
        if K0 % 2:
            j_left, j_right = k - 1, k
            c_u = (
                eberlein(k, j_left, 1)
                + eberlein(k, j_left, k)
                - eberlein(k, j_right, 1)
                - eberlein(k, j_right, k)
            )
            require(c_u == 0, "odd-K endpoint mass cancellation", (K0, c_u))
            value = (
                (eberlein(k, 0, 1) + eberlein(k, 0, k)) * h0
                + (
                    eberlein(k, j_right, 1)
                    + eberlein(k, j_right, k)
                )
                * total
            )
            require(
                value == Fraction(K0 * (4 - K0), s0),
                "odd-K endpoint core identity",
                {"K": K0, "actual": value},
            )
        else:
            j_left, j_right = k - 1, k
            c_u = (
                eberlein(k, j_left, 1)
                - eberlein(k, j_right, 1)
                - 3
            )
            require(c_u == 0, "even-K endpoint mass cancellation", (K0, c_u))
            value = (
                eberlein(k, 0, 1) * h0
                + (eberlein(k, j_right, 1) + 3) * total
            )
            require(
                value == Fraction(K0 * (7 - K0), s0),
                "even-K endpoint core identity",
                {"K": K0, "actual": value},
            )
    print("PASS  core endpoint identity regression checks (K=5..150)")


def main() -> None:
    endpoint_reconstruction()
    interpolation_reconstruction()
    arithmetic_partition()
    sign_domains()
    displayed_independent_family_construction()
    core_endpoint_identities_from_eberlein()
    print("All reconstruction and finite-range checks completed successfully.")


if __name__ == "__main__":
    main()
