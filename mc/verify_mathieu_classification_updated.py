#!/usr/bin/env python3
"""Exact symbolic checks for the revised Mathieu-classification manuscript."""

from __future__ import annotations

import sympy as sp


x, w, z, t, u, y = sp.symbols("x w z t u y")
a, b, c, d = sp.symbols("a b c d")


def constant_term_laurent(expr: sp.Expr, var: sp.Symbol) -> sp.Expr:
    """Return the coefficient of var**0 in an expanded Laurent polynomial."""
    total = sp.Integer(0)
    for term in sp.Add.make_args(sp.expand(expr)):
        power = term.as_powers_dict().get(var, sp.Integer(0))
        if power == 0:
            total += term
    return sp.simplify(total)


# The earlier xz(1,1) witness after the admissible quadratic reparametrization.
f_xz = (1 - w**-1) * ((1 - y) + y * w)
f_xz_weighted = sp.expand(f_xz.subs(y, x**2))
for m in range(1, 7):
    pure = sp.simplify(
        sp.integrate(constant_term_laurent(f_xz_weighted**m, w) * x, (x, 0, 1))
    )
    marked = sp.simplify(
        sp.integrate(
            constant_term_laurent(w**-1 * f_xz_weighted**m, w) * x,
            (x, 0, 1),
        )
    )
    assert pure == 0, (m, pure)
    assert marked == sp.Rational((-1) ** (m - 1), 2 * (m + 1)), (m, marked)


# The explicit abelian pair and its claimed Laurent expansion.
U = 2 * x * (1 - x**2) * w
V = 2 * x * w**-1
T = 1 - 2 * x**2
P_ab = sp.expand((1 + U) * (V - (2 + U) * T**2))
Q_ab = U

claimed_expansion = (
    2 * x * w**-1
    - 2 * (6 * x**4 - 6 * x**2 + 1)
    + 6 * x * (x**2 - 1) * (2 * x**2 - 1) ** 2 * w
    - 4 * x**2 * (x**2 - 1) ** 2 * (2 * x**2 - 1) ** 2 * w**2
)
assert sp.expand(P_ab - claimed_expansion) == 0

spectrum = sorted(
    {
        int(term.as_powers_dict().get(w, sp.Integer(0)))
        for term in sp.Add.make_args(sp.expand(P_ab))
    }
)
assert spectrum == [-1, 0, 1, 2]
assert sp.expand(U * V + T**2 - 1) == 0
assert sp.expand(U * P_ab - (1 + U) * (1 - (1 + U) ** 2 * T**2)) == 0


# Polynomial representatives in the matrix entries and their Mueger--Tuset image.
A0 = a * d - b * c
U0 = -2 * a * b
V0 = 2 * c * d
T0 = a * d + b * c
P0 = (A0 + U0) * (A0**2 * V0 - (2 * A0 + U0) * T0**2)
Q0 = U0

# Invariance under the maximal-torus factor diag(z,z^{-1}).
torus_sub = {a: a * z, b: b * z**-1, c: c * z, d: d * z**-1}
for expr in (A0, U0, V0, T0, P0, Q0):
    assert sp.simplify(expr.subs(torus_sub, simultaneous=True) - expr) == 0

# Square-root-free substitution B(x,w).
B_sub = {
    a: sp.I * w * (1 - x**2),
    b: sp.I * x,
    c: sp.I * x,
    d: -sp.I * w**-1,
}
assert sp.simplify(A0.subs(B_sub, simultaneous=True) - 1) == 0
assert sp.simplify(U0.subs(B_sub, simultaneous=True) - U) == 0
assert sp.simplify(V0.subs(B_sub, simultaneous=True) - V) == 0
assert sp.simplify(T0.subs(B_sub, simultaneous=True) - T) == 0
assert sp.simplify(P0.subs(B_sub, simultaneous=True) - P_ab) == 0
assert sp.simplify(Q0.subs(B_sub, simultaneous=True) - Q_ab) == 0


# Direct exact moment checks for the printed abelian witness.
for m in range(1, 6):
    c_m = sp.integrate((1 - t**2) ** m, (t, 0, 1))
    for s in range(0, m + 3):
        integrand = P_ab**m if s == 0 else Q_ab**s * P_ab**m
        lhs = sp.simplify(2 * sp.integrate(constant_term_laurent(integrand, w) * x, (x, 0, 1)))
        rhs = sp.Integer(0) if s == 0 else sp.simplify(c_m * sp.binomial(m - 1, s - 1))
        assert sp.simplify(lhs - rhs) == 0, (m, s, lhs, rhs)


# Independent reduced Hopf coefficient check through m=8.
for m in range(1, 9):
    c_m = sp.integrate((1 - t**2) ** m, (t, 0, 1))
    for s in range(0, m + 3):
        expr = u ** (s - m) * (1 + u) ** m * (1 - (1 + u) ** 2 * t**2) ** m
        lhs = sp.simplify(sp.integrate(constant_term_laurent(expr, u), (t, 0, 1)))
        rhs = sp.Integer(0) if s == 0 else sp.simplify(c_m * sp.binomial(m - 1, s - 1))
        assert sp.simplify(lhs - rhs) == 0, (m, s, lhs, rhs)

print("All updated exact identities verified, including the weighted xz reparametrization.")
print("Explicit w-spectrum:", spectrum)
