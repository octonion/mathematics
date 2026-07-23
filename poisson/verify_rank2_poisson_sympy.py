#!/usr/bin/env python3
"""Exact SymPy certificate for the proposed PC(2) counterexample.

All arithmetic is over QQ, hence exact.  The script checks:
  * polynomiality and R = x(2-3*x*q),
  * the induced three-variable Jacobian Poisson bracket,
  * the marked-root 3x3 Jacobian,
  * the Hamiltonian-correction vector identity,
  * all six generator Poisson identities,
  * the ordinary 4x4 Jacobian determinant,
  * a polynomial source-coordinate automorphism and its inverse,
  * an explicit three-point collision.
"""

from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# Four-variable construction
# ---------------------------------------------------------------------------
x, q, p, z = sp.symbols("x q p z")
vars4 = (x, q, p, z)

a = 1 - 3*x*q
B = 3*x**2*p + 2*a*z
beta = B - 9*q**2
y = q - x*beta/3
u = x*y

R = 2*x - 3*x**2*y - x**3*beta
S = y + 3*x*(1 + x*y)**2*beta + 3*x*y**2*(4 + 3*x*y)
T = -sp.Rational(1, 2) * (
    (1 + x*y)**3*beta + y**2*(1 + x*y)*(4 + 3*x*y)
)

D0 = sp.Rational(1, 2)*(1 + 3*x*q)*p - 3*q**2*z
H = (
    y**4/sp.Integer(20) * (18*u**2 + 78*u + 125)
    + sp.Rational(3, 10)*beta*y**2*(u**3 + 5*u**2 + 10*u - 5)
    - beta**2/sp.Integer(6)*(9*u + 2)
    - x**2*beta**3/sp.Integer(6)
)
D = D0 + H


def poisson(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return (
        sp.diff(f, p)*sp.diff(g, x) - sp.diff(f, x)*sp.diff(g, p)
        + sp.diff(f, z)*sp.diff(g, q) - sp.diff(f, q)*sp.diff(g, z)
    )


def assert_zero(expr: sp.Expr, generators: tuple[sp.Symbol, ...], label: str) -> None:
    poly = sp.Poly(sp.expand(expr), *generators, domain=sp.QQ)
    assert poly.is_zero, f"FAILED: {label}: {poly.as_expr()}"


def assert_one(expr: sp.Expr, generators: tuple[sp.Symbol, ...], label: str) -> None:
    assert_zero(expr - 1, generators, label)


assert_zero(R - x*(2 - 3*x*q), vars4, "R simplification")

# ---------------------------------------------------------------------------
# Structural three-variable check
# ---------------------------------------------------------------------------
X, Y, W = sp.symbols("X Y W")
vars3 = (X, Y, W)
U = X*Y

r = 2*X - 3*X**2*Y - X**3*W
s = Y + 3*X*(1 + X*Y)**2*W + 3*X*Y**2*(4 + 3*X*Y)
t = -sp.Rational(1, 2) * (
    (1 + X*Y)**3*W + Y**2*(1 + X*Y)*(4 + 3*X*Y)
)
h = (
    Y**4/sp.Integer(20) * (18*U**2 + 78*U + 125)
    + sp.Rational(3, 10)*W*Y**2*(U**3 + 5*U**2 + 10*U - 5)
    - W**2/sp.Integer(6)*(9*U + 2)
    - X**2*W**3/sp.Integer(6)
)

jac_rst = sp.det(sp.Matrix([[sp.diff(f, v) for v in vars3] for f in (r, s, t)]))
assert_zero(jac_rst + 1, vars3, "det d(r,s,t)/d(X,Y,W) = -1")

# Induced brackets on C[x,y,beta]: {f,g} = -Jac(r,f,g).
assert_zero(poisson(x, y) - x**3, vars4, "{x,y}=x^3")
assert_zero(poisson(x, beta) + 3*x**2, vars4, "{x,beta}=-3x^2")
assert_zero(poisson(y, beta) - (-2 + 6*x*y + 3*x**2*beta), vars4,
            "{y,beta}=-2+6xy+3x^2 beta")

# D0 derivation on C[X,Y,W], with q = Y + XW/3.
Q = Y + X*W/3
delta0 = sp.Matrix([
    (1 + 3*X*Q)/2,
    -3*Q**2 - W*(1 + 3*X*Q)/6
             - sp.Rational(3, 2)*X*Q*(W + 21*Q**2),
    sp.Rational(9, 2)*Q*(W + 21*Q**2),
])

grad = lambda f: sp.Matrix([sp.diff(f, v) for v in vars3])
# Exact Hamiltonian-correction identity:
# delta0 + grad(h) x grad(r) + grad(s) x grad(t) = 0.
correction_residual = delta0 + grad(h).cross(grad(r)) + grad(s).cross(grad(t))
for i, component in enumerate(correction_residual):
    assert_zero(component, vars3, f"Hamiltonian correction component {i}")

# ---------------------------------------------------------------------------
# Direct four-variable Poisson and Jacobian certificate
# ---------------------------------------------------------------------------
checks = (
    (poisson(D, R), 1, "{D,R}=1"),
    (poisson(S, T), 1, "{S,T}=1"),
    (poisson(R, S), 0, "{R,S}=0"),
    (poisson(R, T), 0, "{R,T}=0"),
    (poisson(D, S), 0, "{D,S}=0"),
    (poisson(D, T), 0, "{D,T}=0"),
)
for lhs, rhs, label in checks:
    assert_zero(lhs - rhs, vars4, label)

outputs = (R, T, D, S)
J4 = sp.Matrix([[sp.diff(f, v) for v in vars4] for f in outputs])
det4 = J4.det(method="domain-ge")
assert_one(det4, vars4, "det J(R,T,D,S)=1")

# ---------------------------------------------------------------------------
# Source-coordinate automorphism Psi=(x,y,beta,D)
# ---------------------------------------------------------------------------
# Inverse formulas.  E is the new D-coordinate and E0=D0=E-h(X,Y,W).
E = sp.symbols("E")
Qinv = Y + X*W/3
Ainv = 1 - 3*X*Qinv
Cinv = (1 + 3*X*Qinv)/2
E0inv = E - h
Pinv = 3*Qinv**2*(W + 9*Qinv**2) + 2*Ainv*E0inv
Zinv = Cinv*(W + 9*Qinv**2) - 3*X**2*E0inv

forward_sub = {X: x, Y: y, W: beta, E: D}
assert_zero(Qinv.subs(forward_sub) - q, vars4, "Psi inverse: q")
assert_zero(Pinv.subs(forward_sub) - p, vars4, "Psi inverse: p")
assert_zero(Zinv.subs(forward_sub) - z, vars4, "Psi inverse: z")

Jpsi = sp.Matrix([[sp.diff(f, v) for v in vars4] for f in (x, y, beta, D)])
assert_zero(Jpsi.det(method="domain-ge") + 1, vars4, "det J Psi=-1")

# ---------------------------------------------------------------------------
# Explicit three-point collision (coordinates are (x,q,p,z))
# ---------------------------------------------------------------------------
points = (
    (sp.Rational(0), sp.Rational(0), sp.Rational(1, 24), sp.Rational(-1, 8)),
    (sp.Rational(1), sp.Rational(2, 3), sp.Rational(247, 96), sp.Rational(-89, 64)),
    (sp.Rational(-1), sp.Rational(-2, 3), sp.Rational(247, 96), sp.Rational(-89, 64)),
)
target = (sp.Rational(0), sp.Rational(1, 8), sp.Rational(0), sp.Rational(0))
for i, point in enumerate(points, start=1):
    image = tuple(sp.cancel(f.subs(dict(zip(vars4, point)))) for f in outputs)
    assert image == target, f"FAILED: collision point {i} maps to {image}"

# Polynomial size data are useful for detecting transcription mistakes.
size_data = {}
for name, expr in (("beta", beta), ("y", y), ("R", R), ("S", S),
                   ("T", T), ("H", H), ("D", D)):
    poly = sp.Poly(sp.expand(expr), *vars4, domain=sp.QQ)
    size_data[name] = (len(poly.terms()), poly.total_degree())

print(f"SymPy {sp.__version__}")
print("All exact assertions passed.")
print("Term counts and total degrees:", size_data)
print("det d(r,s,t)/d(X,Y,W) =", sp.expand(jac_rst))
print("det J(R,T,D,S) =", sp.expand(det4))
print("Three collision points map to", target)
