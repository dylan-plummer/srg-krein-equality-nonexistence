#!/usr/bin/env python3
"""Read the machine-readable BCN tables (IA3P, IA4a, IA5a, IA6b,...,IA17b).
Usage: python distance_arrays.py PATH_TO_DIRECTORY
All rational restricted eigenspaces are tested, including those of imprimitive graphs.
"""
import re
import sys
from collections import Counter
from fractions import Fraction as Q
from math import gcd, lcm
from pathlib import Path

from frame_test import obstruction

CATALOG = ["IA3P", "IA4a", "IA5a", *(f"IA{i}b" for i in range(6, 18))]


def roots(b, c):
    k = b[0]
    older, p = [], [1]
    for i in range(len(b)):
        a = k-b[i]-c[i]
        q = [0, *p]
        for j, coeff in enumerate(p):
            q[j] -= a*coeff
        for j, coeff in enumerate(older):
            q[j] -= b[i-1]*c[i]*coeff
        older, p = p, q
    # Rational eigenvalues of an integral adjacency matrix are integers.
    while p[0] == 0:
        yield 0
        p = p[1:]
    for x in range(-k, k):   # excludes the trivial eigenvalue k
        if not x or p[0] % x:
            continue
        value = 0
        for coeff in reversed(p):
            value = value*x+coeff
        if not value:
            yield x


def frames(up, down):
    b, c = [*up, 0], [0, *down]
    k = b[0]
    counts = valencies(up, down)
    n = sum(counts)
    for theta in roots(b, c):
        t = [Q(1)]
        for i in range(len(up)):
            a = k-b[i]-c[i]
            t.append(((theta-a)*t[i]-(c[i]*t[i-1] if i else 0))/b[i])
        assert sum(h*z for h,z in zip(counts,t)) == 0
        cubic = sum(h*z**3 for h,z in zip(counts,t)) == 0
        g0 = lcm(*(v.denominator for v in t))
        g = [int(v*g0) for v in t]
        factor = sum(h*x*x for h,x in zip(counts,g)) / Q(g0)
        dimension = n*g0/factor
        assert factor.denominator == dimension.denominator == 1 and gcd(*g) == 1
        yield n, theta, int(dimension), int(factor), g, cubic


def valencies(up, down):
    counts = [1]
    for b, c in zip(up, down):
        val, rem = divmod(counts[-1]*b, c)
        assert rem == 0
        counts.append(val)
    return counts


def controls():
    # The indicated distance fusions must both be the triangle-free (4732,166,0,6) SRG.
    for up, down, distances in [((165,160,6,1), (1,6,160,165), (3,4)),
                                ((415,390,16), (1,39,400), (3,))]:
        k = valencies(up, down)
        assert sum(k[i] for i in distances) == 166
        spectrum = Counter()
        for n, theta, m, frame, g, cubic in frames(up, down):
            ev = Q(sum(k[i]*g[i] for i in distances), g[0])
            spectrum[ev] += m
        assert spectrum == {10: 2905, -16: 1826}
    # Actual Hamming graphs H(d,q) and Johnson graphs J(v,d).
    for d in range(2, 8):
        for q in range(2, 9):
            hamming = [(d-i)*(q-1) for i in range(d)], list(range(1,d+1))
            v = 2*d+q-2
            johnson = [(d-i)*(v-d-i) for i in range(d)], [i*i for i in range(1,d+1)]
            for arr in (hamming, johnson):
                records = list(frames(*arr))
                assert len(records) == d
                for n, theta, m, frame, g, cubic in records:
                    assert not obstruction(n, m, frame, g, cubic), (arr, theta)


def main(directory):
    seen = set()
    excluded = set()
    counts = Counter()
    controls()
    print("table", "intersection_array", "n", "theta", "m", "frame", "gram", "zero_cubic", "test", sep="\t")
    for name in CATALOG:
        raw = (directory/name).read_text()
        for a, c in re.findall(r'i\(([\d,\s]+);([\d,\s]+)\)', raw):
            arr = tuple(map(int,a.split(','))), tuple(map(int,c.split(',')))
            if arr in seen:
                continue
            seen.add(arr)
            for n, theta, m, frame, g, cubic in frames(*arr):
                counts["spaces"] += 1
                counts["cubics"] += cubic
                reason = obstruction(n,m,frame,g,cubic)
                if reason:
                    excluded.add(arr)
                    print(name, arr, n, theta, m, frame, g, cubic, reason, sep='\t')
                    counts["exclusions"] += 1
    print("arrays", len(seen), "rational eigenspaces", counts["spaces"], "rational cubic-zero eigenspaces", counts["cubics"],
          "excluded eigenspaces", counts["exclusions"], "excluded arrays", len(excluded), file=sys.stderr)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
