#!/usr/bin/env python3
"""Independent exact sparse-polynomial audit using only Python's standard library.

This deliberately does not import SymPy or any other CAS.  Coefficients are
fractions.Fraction values and polynomials are dictionaries of monomials.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from typing import Dict, Iterable, Tuple, Union

Monomial = Tuple[int, int, int, int]
Scalar = Union[int, Fraction]


class Poly4:
    __slots__ = ("terms",)

    def __init__(self, terms: Dict[Monomial, Fraction] | None = None) -> None:
        self.terms = {m: Fraction(c) for m, c in (terms or {}).items() if c}

    @staticmethod
    def const(c: Scalar) -> "Poly4":
        c = Fraction(c)
        return Poly4({(0, 0, 0, 0): c}) if c else Poly4()

    @staticmethod
    def var(index: int) -> "Poly4":
        exponents = [0, 0, 0, 0]
        exponents[index] = 1
        return Poly4({tuple(exponents): Fraction(1)})

    def __add__(self, other: object) -> "Poly4":
        rhs = as_poly(other)
        result = dict(self.terms)
        for monomial, coefficient in rhs.terms.items():
            result[monomial] = result.get(monomial, Fraction(0)) + coefficient
            if not result[monomial]:
                del result[monomial]
        return Poly4(result)

    __radd__ = __add__

    def __neg__(self) -> "Poly4":
        return Poly4({m: -c for m, c in self.terms.items()})

    def __sub__(self, other: object) -> "Poly4":
        return self + (-as_poly(other))

    def __rsub__(self, other: object) -> "Poly4":
        return as_poly(other) - self

    def __mul__(self, other: object) -> "Poly4":
        rhs = as_poly(other)
        if not self.terms or not rhs.terms:
            return Poly4()
        result: Dict[Monomial, Fraction] = {}
        for m, c in self.terms.items():
            for n, d in rhs.terms.items():
                monomial = tuple(m[i] + n[i] for i in range(4))
                result[monomial] = result.get(monomial, Fraction(0)) + c*d
        return Poly4({m: c for m, c in result.items() if c})

    __rmul__ = __mul__

    def __truediv__(self, scalar: Scalar) -> "Poly4":
        scalar = Fraction(scalar)
        if not scalar:
            raise ZeroDivisionError
        return Poly4({m: c/scalar for m, c in self.terms.items()})

    def __pow__(self, exponent: int) -> "Poly4":
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("Exponent must be a nonnegative integer")
        result = Poly4.const(1)
        base = self
        n = exponent
        while n:
            if n & 1:
                result = result*base
            n >>= 1
            if n:
                base = base*base
        return result

    def diff(self, index: int) -> "Poly4":
        result: Dict[Monomial, Fraction] = {}
        for monomial, coefficient in self.terms.items():
            exponent = monomial[index]
            if exponent:
                new_monomial = list(monomial)
                new_monomial[index] -= 1
                result[tuple(new_monomial)] = coefficient*exponent
        return Poly4(result)

    def evaluate(self, point: Iterable[Scalar]) -> Fraction:
        values = tuple(Fraction(v) for v in point)
        if len(values) != 4:
            raise ValueError("Expected four coordinates")
        total = Fraction(0)
        for monomial, coefficient in self.terms.items():
            term = coefficient
            for value, exponent in zip(values, monomial):
                term *= value**exponent
            total += term
        return total

    def total_degree(self) -> int:
        return max((sum(m) for m in self.terms), default=-1)

    def __eq__(self, other: object) -> bool:
        try:
            rhs = as_poly(other)
        except TypeError:
            return False
        return self.terms == rhs.terms


def as_poly(value: object) -> Poly4:
    if isinstance(value, Poly4):
        return value
    if isinstance(value, (int, Fraction)):
        return Poly4.const(value)
    raise TypeError(f"Cannot convert {type(value)!r} to Poly4")


def poisson(f: Poly4, g: Poly4) -> Poly4:
    # Variable order: (x,q,p,z).
    return (
        f.diff(2)*g.diff(0) - f.diff(0)*g.diff(2)
        + f.diff(3)*g.diff(1) - f.diff(1)*g.diff(3)
    )


def permutation_sign(perm: Tuple[int, ...]) -> int:
    inversions = sum(
        perm[i] > perm[j]
        for i in range(len(perm))
        for j in range(i + 1, len(perm))
    )
    return -1 if inversions % 2 else 1


def determinant(matrix: Tuple[Tuple[Poly4, ...], ...]) -> Poly4:
    n = len(matrix)
    result = Poly4.const(0)
    for perm in permutations(range(n)):
        term = Poly4.const(permutation_sign(perm))
        for row, column in enumerate(perm):
            term = term*matrix[row][column]
        result = result + term
    return result


x, q, p, z = (Poly4.var(i) for i in range(4))
a = 1 - 3*x*q
B = 3*x**2*p + 2*a*z
beta = B - 9*q**2
y = q - x*beta/3
u = x*y

R = 2*x - 3*x**2*y - x**3*beta
S = y + 3*x*(1 + x*y)**2*beta + 3*x*y**2*(4 + 3*x*y)
T = -((1 + x*y)**3*beta + y**2*(1 + x*y)*(4 + 3*x*y))/2
D0 = (1 + 3*x*q)*p/2 - 3*q**2*z
H = (
    y**4*(18*u**2 + 78*u + 125)/20
    + 3*beta*y**2*(u**3 + 5*u**2 + 10*u - 5)/10
    - beta**2*(9*u + 2)/6
    - x**2*beta**3/6
)
D = D0 + H

assert R == x*(2 - 3*x*q)

checks = (
    (poisson(D, R), 1, "{D,R}=1"),
    (poisson(S, T), 1, "{S,T}=1"),
    (poisson(R, S), 0, "{R,S}=0"),
    (poisson(R, T), 0, "{R,T}=0"),
    (poisson(D, S), 0, "{D,S}=0"),
    (poisson(D, T), 0, "{D,T}=0"),
)
for actual, expected, label in checks:
    assert actual == expected, f"FAILED: {label}; {len(actual.terms)} residual terms"

outputs = (R, T, D, S)
jacobian = tuple(tuple(f.diff(j) for j in range(4)) for f in outputs)
det_jacobian = determinant(jacobian)
assert det_jacobian == 1

points = (
    (Fraction(0), Fraction(0), Fraction(1, 24), Fraction(-1, 8)),
    (Fraction(1), Fraction(2, 3), Fraction(247, 96), Fraction(-89, 64)),
    (Fraction(-1), Fraction(-2, 3), Fraction(247, 96), Fraction(-89, 64)),
)
target = (Fraction(0), Fraction(1, 8), Fraction(0), Fraction(0))
for index, point in enumerate(points, start=1):
    image = tuple(f.evaluate(point) for f in outputs)
    assert image == target, f"FAILED: point {index} maps to {image}"

size_data = {
    name: (len(poly.terms), poly.total_degree())
    for name, poly in (
        ("beta", beta), ("y", y), ("R", R), ("S", S),
        ("T", T), ("H", H), ("D", D),
    )
}

print("All independent sparse-polynomial assertions passed.")
print("Term counts and total degrees:", size_data)
print("det J(R,T,D,S) =", det_jacobian.terms)
print("Three collision points map to", target)
