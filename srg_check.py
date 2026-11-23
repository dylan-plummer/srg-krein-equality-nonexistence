#!/usr/bin/env python3
"""Independent check of the infinite families of Section 7 of the paper,
using the generic criterion of frame_check.py (Theorem 3.1 + Lemma 3.2 of the paper)
applied to b(rI-A+uJ) in both orientations of each parameter set.

Usage: python -B srg_check.py   (run from this directory)
"""
import sys
from fractions import Fraction as Fr
from math import gcd, isqrt

import sympy as sp

from frame_check import criterion


def orient(n, k, lam, mu):
    """Return the frame data (n, m, F, g, cubic) for the s-eigenspace, or None."""
    D = (lam - mu) ** 2 + 4 * (k - mu)
    t = isqrt(D)
    if t * t != D or (lam - mu + t) % 2:
        return None
    r, s = (lam - mu + t) // 2, (lam - mu - t) // 2
    E = r - s
    u = Fr(k - r, n)
    a, b = u.numerator, u.denominator
    m = Fr(n * (b * r + a), b * E)
    if m.denominator != 1:
        return None
    cubic = E * (1 - 2 * u) == r * (r + 1)
    return n, int(m), b * E, [b * r + a, a, a - b], cubic


def excluded(row):
    n, k, lam, mu = row
    comp = (n, n - k - 1, n - 2 - 2 * k + mu, n - 2 * k + lam)
    out = []
    for R in (row, comp):
        o = orient(*R)
        if o:
            v = criterion(*o)
            if v:
                out.append((R, v))
    return out


def ok_primes(q):
    return all(p % 4 == 1 for p in sp.factorint(q))


def basic_feasible(n, k, lam, mu):
    if k * (k - lam - 1) != (n - k - 1) * mu:
        return False
    o = orient(n, k, lam, mu)
    return o is not None


def main():
    # (A) Pell family
    x, q = 1, 1
    pell = 0
    for _ in range(24):
        x, q = 3 * x + 4 * q, 2 * x + 3 * q
        if q % 8 != 5:
            continue
        c = (1 + x) // 2 if (1 + x) // 2 % 2 else (1 - x) // 2
        row = (2 * q * q, q * q - c, c * (c - 3), c * (c - 1))
        assert basic_feasible(*row), row
        assert excluded(row), row
        pell += 1
    print("Pell family: %d members (q==5 mod 8) excluded, first q=%d" % (pell, 5), file=sys.stderr)

    # (B) fixed (lambda, mu) families n = q^2(q^2-1)/mu, k = q^2-c
    classes = {-1: (8, {5}), 3: (24, {13, 17}), -2: (24, {13, 17}), 4: (48, {17}), -3: (48, {17})}
    cnt = 0
    for c, (mod, res) in classes.items():
        mu = c * (c - 1)
        for q in range(5, 3000, 2):
            if q % mod not in res or not ok_primes(q):
                continue
            N = q * q * (q * q - 1)
            if N % mu:
                continue
            row = (N // mu, q * q - c, c * (c - 3), mu)
            if not basic_feasible(*row):
                continue
            assert excluded(row), (c, q, row)
            cnt += 1
    print("fixed (lambda,mu) families: %d instances q<3000, all excluded" % cnt, file=sys.stderr)

    # (C) triangle-free mu=2 (k=q^2+1) and mu=4 (k=q^2), q==5 mod 8, primes 1 mod 4
    cnt = 0
    for q in range(5, 3000, 8):
        if not ok_primes(q):
            continue
        for mu, k in ((2, q * q + 1), (4, q * q)):
            if (k * (k - 1)) % mu:
                continue
            row = (1 + k + k * (k - 1) // mu, k, 0, mu)
            if not basic_feasible(*row):
                continue
            assert excluded(row), row
            cnt += 1
    print("triangle-free mu=2,4: %d instances q<3000, all excluded" % cnt, file=sys.stderr)

    # (D) named rows
    for row in [(1682, 820, 378, 420), (4732, 166, 0, 6), (4732, 171, 10, 6), (14196, 170, 4, 2),
                (13872, 286, 0, 6), (13872, 291, 10, 6), (6936, 285, 4, 12), (6936, 292, 18, 12),
                (14536, 170, 0, 2), (7268, 169, 0, 4), (352, 26, 0, 2), (176, 25, 0, 4),
                (24167, 5590, 1053, 1365), (181545, 25960, 2975, 3835), (421685, 76704, 12478, 14280),
                (595051, 245050, 98049, 102921), (731445, 332500, 147275, 154375),
                (1691285, 178080, 15950, 19080), (4001885, 1740000, 745550, 765000)]:
        e = excluded(row)
        print(row, "EXCLUDED" if e else "passes", e[0][1] if e else "", sep="\t")


if __name__ == "__main__":
    main()
