"""Exact congruences for a rational constant-diagonal tight frame."""
from fractions import Fraction
from math import gcd


def primes(n):
    result = []
    p = 2
    while p*p <= n:
        if n % p == 0:
            result.append(p)
            while n % p == 0:
                n //= p
        p = 3 if p == 2 else p+2
    return result + ([n] if n > 1 else [])


def v2(n):
    return (n & -n).bit_length()-1


def obstruction(n, m, frame, values, cubic=False):
    """Gram values are primitive integral, diagonal first, and the frame sums to zero."""
    odd = values[0] % 2
    modulus = 16 if cubic and odd else 8
    target = n if odd else 0
    if target % gcd(frame, modulus):
        return "integrality"
    if v2(frame) > v2(modulus)-2:
        return None
    q = Fraction(target, frame)
    norm = q.numerator * pow(q.denominator, -1, 4) % 4
    if (norm-m) % 2:
        return "dimension"
    b = gcd(*(abs(x-values[0]) for x in values))
    unknown = primes(Fraction(frame, b*b if cubic else b).numerator)
    if 2 not in unknown and (m+1)*v2(b) % 2:
        return "determinant valuation"
    if all(p % 4 != 3 for p in unknown):
        disc_unit = pow(b >> v2(b), m+1, 4)
        if (m+1-disc_unit-norm) % 4:
            return "determinant unit"
    return None
