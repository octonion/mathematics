import sympy as sp

x, q, p, z = sp.symbols('x q p z')

# Canonical bracket: {p,x}=1 and {z,q}=1.
def PB(f, g):
    return sp.expand(
        sp.diff(f, p) * sp.diff(g, x)
        - sp.diff(f, x) * sp.diff(g, p)
        + sp.diff(f, z) * sp.diff(g, q)
        - sp.diff(f, q) * sp.diff(g, z)
    )

a = 1 - 3*x*q
B = 3*x**2*p + 2*a*z
beta = B - 9*q**2
y = q - x*beta/3

R = 2*x - 3*x**2*y - x**3*beta
S = y + 3*x*(1 + x*y)**2*beta + 3*x*y**2*(4 + 3*x*y)
T = -sp.Rational(1, 2) * (
    (1 + x*y)**3*beta + y**2*(1 + x*y)*(4 + 3*x*y)
)

u = x*y
H = (
    sp.Rational(1, 20)*y**4*(18*u**2 + 78*u + 125)
    + sp.Rational(3, 10)*beta*y**2*(u**3 + 5*u**2 + 10*u - 5)
    - sp.Rational(1, 6)*beta**2*(9*u + 2)
    - sp.Rational(1, 6)*x**2*beta**3
)
D0 = sp.Rational(1, 2)*(1 + 3*x*q)*p - 3*q**2*z
D = sp.expand(D0 + H)

identities = {
    '{D,R}-1': PB(D, R) - 1,
    '{S,T}-1': PB(S, T) - 1,
    '{R,T}': PB(R, T),
    '{R,S}': PB(R, S),
    '{D,T}': PB(D, T),
    '{D,S}': PB(D, S),
}

for name, expression in identities.items():
    value = sp.factor(expression)
    print(f'{name}: {value}')
    assert value == 0

# Generator order is (x,q,p,z) -> (R,T,D,S).
J = sp.Matrix([
    [sp.diff(F, v) for v in (x, q, p, z)]
    for F in (R, T, D, S)
])
print('Jacobian determinant:', sp.factor(J.det()))
assert sp.factor(J.det()) == 1

print('R factorization:', sp.factor(R))
assert sp.expand(R - x*(2 - 3*x*q)) == 0
print('All exact checks passed.')
