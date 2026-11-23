#!/usr/bin/env python3
"""Independent implementation of the frame criterion (Theorem 3.1 + Lemma 3.2
of the paper) for eigenspaces of distance-regular graphs, written without
reference to frame_test.py.

Usage:
  python -I frame_check.py               # controls on known DRG families
  python -I frame_check.py BCN_DIR       # also scan the BCN array files
"""
import re
import sys
from fractions import Fraction as Fr
from math import gcd
from pathlib import Path

import numpy as np
import sympy as sp


def v2(x):
    x = abs(x)
    return (x & -x).bit_length() - 1


def odd_primes(n):
    return [p for p in sp.factorint(abs(n)) if p != 2]


def eigenspaces(b, c):
    """b=[b0..b_{d-1}], c=[c1..cd]; yields (n,theta,m,F,g,cubic) for rational theta != k."""
    d = len(b)
    B = list(b) + [0]
    C = [0] + list(c)
    k = b[0]
    kk = [1]
    for i in range(d):
        q, r = divmod(kk[-1] * B[i], C[i + 1])
        if r:
            return
        kk.append(q)
    n = sum(kk)
    T = np.zeros((d + 1, d + 1))
    for i in range(d + 1):
        T[i, i] = k - B[i] - C[i]
        if i < d:
            T[i, i + 1] = B[i]
        if i > 0:
            T[i, i - 1] = C[i]
    # Rational eigenvalues of an integral matrix are integers; T has real
    # spectrum (it is similar to a symmetric matrix).  Round, then verify exactly.
    cands = sorted({int(round(z.real)) for z in np.linalg.eigvals(T)})
    for theta in cands:
        if theta == k:
            continue
        u = [Fr(1), Fr(theta, k)]
        for i in range(1, d):
            a_i = k - B[i] - C[i]
            u.append(((theta - a_i) * u[i] - C[i] * u[i - 1]) / B[i])
        u = u[:d + 1]
        # theta is an eigenvalue iff the last row of the recurrence holds
        if C[d] * u[d - 1] + (k - C[d]) * u[d] != theta * u[d]:
            continue
        den = 1
        for z in u:
            den = den * z.denominator // gcd(den, z.denominator)
        g = [int(z * den) for z in u]
        h = 0
        for z in g:
            h = gcd(h, z)
        g = [z // h for z in g]
        m = Fr(n * g[0] ** 2, sum(kh * gh * gh for kh, gh in zip(kk, g)))
        if m.denominator != 1:
            continue
        m = int(m)
        F = Fr(g[0] * n, m)
        assert F.denominator == 1
        cubic = sum(kh * gh ** 3 for kh, gh in zip(kk, g)) == 0
        yield n, theta, m, int(F), g, cubic


def criterion(n, m, F, g, cubic):
    """Return None (no contradiction / not applicable) or a reason string."""
    g0 = g[0]
    if g0 % 2:
        T, M = n, (16 if cubic else 8)
    else:
        T, M = 0, 8          # first moment vanishes since theta != k
    if v2(F) > v2(M) - 2:
        return None
    # F w^2 == T (mod M) with w^2 in Z_(2): needs v2(T) >= v2(F)
    if T and v2(T) < v2(F):
        return "integrality"
    t = Fr(T, F)
    t4 = (t.numerator * pow(t.denominator, -1, 4)) % 4 if T else 0
    c = 0
    for z in g:
        c = gcd(c, z - g0)
    c = abs(c)
    p = Fr(F, c * c if cubic else c).numerator
    # Determinant square class = c^(m+1) * (primes dividing p), rationally.
    if p % 2 and ((m + 1) * v2(c)) % 2:
        return "det valuation odd"
    oc = c >> v2(c)
    D = pow(oc, m + 1, 4)
    if (m + 1 - D - t4) % 4 == 0:
        return None
    # The congruence fails if D is as computed; D is determined only when
    # every odd prime dividing p is 1 mod 4 (factor only in this rare case).
    if any(q % 4 == 3 for q in odd_primes(p)):
        return None
    return "char norm %d vs m+1-D=%d" % (t4, (m + 1 - D) % 4)


def q_int(q, j):
    return (q ** j - 1) // (q - 1)


def controls():
    arrays = []
    for d in range(2, 7):
        for q in range(2, 8):
            arrays.append(([(d - i) * (q - 1) for i in range(d)], list(range(1, d + 1))))
        for v in range(2 * d, 2 * d + 9):
            arrays.append(([(d - i) * (v - d - i) for i in range(d)], [i * i for i in range(1, d + 1)]))
    for q in (2, 3, 4, 5):
        for d in range(2, 5):
            for e in (0, 1, 2):                                   # dual polar
                arrays.append(([q ** (i + e) * q_int(q, d - i) for i in range(d)],
                               [q_int(q, i) for i in range(1, d + 1)]))
            for nn in range(2 * d, 2 * d + 3):                   # Grassmann
                arrays.append(([q ** (2 * i + 1) * q_int(q, d - i) * q_int(q, nn - d - i) for i in range(d)],
                               [q_int(q, i) ** 2 for i in range(1, d + 1)]))
            for e in range(d, d + 3):                             # bilinear forms
                arrays.append(([(q ** d - q ** i) * (q ** e - q ** i) // (q - 1) for i in range(d)],
                               [q ** (i - 1) * q_int(q, i) for i in range(1, d + 1)]))
        for (s, t) in ((q, q), (q, q ** 3), (q ** 3, q)):          # generalized hexagons
            arrays.append(([s * (t + 1), s * t, s * t], [1, 1, t + 1]))
    for q in (2, 8, 32):                                          # generalized octagons
        for (s, t) in ((q, q * q), (q * q, q)):
            arrays.append(([s * (t + 1), s * t, s * t, s * t], [1, 1, 1, t + 1]))
    for nn in range(4, 16):                                       # halved cubes
        dd = nn // 2
        arrays.append(([(nn - 2 * i) * (nn - 2 * i - 1) // 2 for i in range(dd)],
                       [i * (2 * i - 1) for i in range(1, dd + 1)]))
    arrays += [  # sporadic known graphs
        ([3, 2, 2, 2, 1, 1, 1], [1, 1, 1, 1, 1, 1, 3]),           # Biggs-Smith
        ([3, 2, 2, 2, 2, 1, 1, 1], [1, 1, 1, 1, 2, 2, 2, 3]),     # Foster
        ([3, 2, 2, 1], [1, 1, 1, 2]),                             # Coxeter
        ([5, 4, 1, 1], [1, 1, 4, 5]),                             # Wells
        ([10, 6, 4], [1, 2, 5]),                                  # Doro
        ([10, 8, 8, 2], [1, 1, 4, 5]),                            # Hall-Janko near octagon
        ([11, 10, 6, 1], [1, 1, 5, 11]),                          # Livingstone
        ([30, 28, 24], [1, 3, 15]),                               # M24 on 759
        ([7, 6, 4, 4], [1, 1, 1, 6]),                             # doubly truncated Witt graph, 330 vertices
        ([22, 21, 20], [1, 2, 6]),                                # coset graph of truncated binary Golay code, 1024 vertices
    ]
    tested = spaces = 0
    for b, c in arrays:
        for n, theta, m, F, g, cubic in eigenspaces(b, c):
            spaces += 1
            r = criterion(n, m, F, g, cubic)
            assert r is None, (b, c, theta, r)
        tested += 1
    print("controls: %d arrays, %d rational eigenspaces, 0 contradictions" % (tested, spaces), file=sys.stderr)


def scan(directory):
    seen = set()
    hits = []
    for f in sorted(Path(directory).glob("IA*")):
        if f.suffix:
            continue
        for a, c in re.findall(r'i\(([\d,\s]+);([\d,\s]+)\)', f.read_text()):
            arr = (tuple(map(int, a.split(','))), tuple(map(int, c.split(','))))
            if arr in seen:
                continue
            seen.add(arr)
            for n, theta, m, F, g, cubic in eigenspaces(*arr):
                r = criterion(n, m, F, g, cubic)
                if r:
                    hits.append((n, arr, theta, m, F, g, cubic, r))
    for h in sorted(hits, key=lambda z: (z[0], z[1])):
        print(*h, sep="\t")
    print("arrays %d, failing eigenspaces %d, failing arrays %d" %
          (len(seen), len(hits), len({h[1] for h in hits})), file=sys.stderr)


if __name__ == "__main__":
    controls()
    if len(sys.argv) > 1:
        scan(sys.argv[1])
