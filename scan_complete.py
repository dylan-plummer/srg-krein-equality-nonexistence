#!/usr/bin/env python3
"""Complete enumeration of Smith parameter sets (q^s_ss = 0, integral eigenvalues)
with positive eigenvalue r <= RMAX, and every test of the paper.

For r >= 2, with P = r(r+1) and E = P + g, integrality of n forces
g | 4 r^3 (r+1)^3 (r-1);  for r = 1, (E+2) | 24.  So each r gives finitely many sets.

Columns: n k lambda mu r s u m feasible checks A V S2P S2S tight4 realized
  feasible = parent passes integrality, Krein and absolute bounds
  checks   = "passes the stated feasibility checks": feasible, integral forced local gaps,
             and each nontrivial subconstituent feasible or a disjoint union of equal
             cliques (or the complement of one)
  A   = Theorem 1.1 verdict (excl / pass / undet / na), using the neighbor basis when lambda = 0
  V   = Theorem 6.1 verdict
  S2P = second-moment test on the graph (both orientations) excludes it
  S2S = second-moment test on a subconstituent (both orientations) excludes it
Usage: python -B scan_complete.py RMAX OUT.tsv      (run from this directory)
"""
import csv
import sys
from pathlib import Path
from fractions import Fraction as Fr
from math import gcd, isqrt

import sympy as sp

from frame_check import criterion

RMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 135
OUT = sys.argv[2] if len(sys.argv) > 2 else None


def v2(x):
    x = abs(x)
    return (x & -x).bit_length() - 1


def oddpart(x):
    x = Fr(x)
    a, b = abs(x.numerator), x.denominator
    return (a >> v2(a)) * pow(b >> v2(b), -1, 8)


def bad_primes(x):
    return [p for p in sp.factorint(abs(x)) if p % 4 == 3]


def eig(n, k, lam, mu):
    D = (lam - mu) ** 2 + 4 * (k - mu)
    t = isqrt(D)
    if t * t != D or (lam - mu + t) % 2:
        return None
    r, s = (lam - mu + t) // 2, (lam - mu - t) // 2
    f = Fr(-k - (n - 1) * s, r - s)   # multiplicity of r
    g = n - 1 - f                     # multiplicity of s
    if f.denominator != 1 or f <= 0 or g <= 0:
        return None
    return r, s, int(f), int(g)


def feasible(n, k, lam, mu):
    if not (0 < mu < k < n - 1 and lam >= 0) or k * (k - lam - 1) != (n - k - 1) * mu:
        return False
    e = eig(n, k, lam, mu)
    if e is None:
        return False
    r, s, f, g = e
    if (s + 1) * (k + s + 2 * r * s) > (k + s) * (r + 1) ** 2:
        return False
    if (r + 1) * (k + r + 2 * r * s) > (k + r) * (s + 1) ** 2:
        return False
    if 2 * n > f * (f + 3) or 2 * n > g * (g + 3):
        return False
    return True


def local_ok(order, deg, cl, cm):
    """subconstituent check: feasible SRG, or a disjoint union of equal cliques / its complement."""
    if cm == 0:                          # union of cliques K_{deg+1}
        return cl == deg - 1 and order % (deg + 1) == 0
    if cm == deg:                        # complete multipartite: complement is a union of cliques
        return local_ok(order, order - deg - 1, order - 2 - 2 * deg + cm, 0)
    return feasible(order, deg, cl, cm)


def passes_checks(n, k, lam, mu, r, s):
    """the single feasibility predicate used by the census, the tables and the plot."""
    if not feasible(n, k, lam, mu):
        return False
    E, u = r - s, Fr(k - r, n)
    for order, deg, gap in ((k, lam, E * u), (n - k - 1, k - mu, E * (1 - u))):
        if deg not in (0, order - 1) and gap.denominator != 1:
            return False
    return all(c is None or local_ok(*c) for c in subconstituents(n, k, lam, mu, r, s))


def frame_s(n, k, lam, mu):
    """frame data (n, m, F, gram, cubic) for the s-eigenspace of b(rI-A+uJ)."""
    e = eig(n, k, lam, mu)
    if e is None:
        return None
    r, s, f, m = e
    E = r - s
    u = Fr(k - r, n)
    a, b = u.numerator, u.denominator
    cubic = E * (1 - 2 * u) == r * (r + 1)
    return n, m, b * E, [b * r + a, a, a - b], cubic


def comp(row):
    n, k, lam, mu = row
    return (n, n - k - 1, n - 2 - 2 * k + mu, n - 2 * k + lam)


def second_moment(row):
    for R in (row, comp(row)):
        fr = frame_s(*R)
        if fr and criterion(fr[0], fr[1], fr[2], fr[3], False):
            return True
    return False


def subconstituents(n, k, lam, mu, r, s):
    E = r - s
    u = Fr(k - r, n)
    out = []
    for order, deg, gap in ((k, lam, E * u), (n - k - 1, k - mu, E * (1 - u))):
        if deg == 0 or deg == order - 1 or gap.denominator != 1:
            out.append(None)
            continue
        ls = r - int(gap)
        cm = deg + r * ls
        cl = cm + r + ls
        out.append((order, deg, cl, cm))
    return out


def thm_main(n, k, lam, mu, r, s, m):
    E = r - s
    u = Fr(k - r, n)
    a, b = u.numerator, u.denominator
    h = b * r + a
    F = b * E
    if lam == 0 and m == k:          # Kre(r): neighbor basis, Gram brI + aJ
        if h % 2 == 0 or v2(F) > 2:
            return 'na'
        det = Fr(b * r) ** ((k - 1) % 2) * (b * r + a * k)
        D = oddpart(det) % 4
        t = Fr(n, F)
        t4 = (t.numerator * pow(t.denominator, -1, 4)) % 4
        return 'excl' if (t4 - (m + 1 - D)) % 4 else 'pass'
    if (h % 2 and v2(F) > 2) or (h % 2 == 0 and v2(F) > 1):
        return 'na'
    if criterion(n, m, F, [h, a, a - b], True):
        return 'excl'
    L = Fr(E, b).numerator
    return 'undet' if any(p != 2 for p in bad_primes(L)) else 'pass'


def thm_vertex(n, k, lam, mu, r, s, m):
    E = r - s
    u = Fr(k - r, n)
    a, b = u.numerator, u.denominator
    if b % 2 or E % b:
        return 'na'
    L = E // b
    if L % 4 == 0:
        return 'na'
    if m % 2 == 0:
        return 'excl'
    if bad_primes(L):
        return 'undet'
    return 'excl' if (b * r + a) % 4 != (-1) ** ((m - 1) // 2) % 4 else 'pass'


def smith_sets(rmax):
    for r in range(1, rmax + 1):
        P = r * (r + 1)
        if r == 1:
            gs = [E - P for E in range(P + 1, 23) if 24 % (E + 2) == 0]
        else:
            gs = sp.divisors(4 * r ** 3 * (r + 1) ** 3 * (r - 1))
        for g in gs:
            E = P + g
            num = 2 * E * E * (E * (2 * r + 1) - 3 * P)
            den = g * (2 * P + g)
            if num % den:
                continue
            n = num // den
            kf = r + Fr(n * g, 2 * E)
            if kf.denominator != 1:
                continue
            k = int(kf)
            s = r - E
            mu = k + r * s
            lam = mu + r + s
            if not (0 < mu < k < n - 1 and lam >= 0):
                continue
            e = eig(n, k, lam, mu)
            if e is None or e[0] != r:
                continue
            yield n, k, lam, mu, r, s, e[3]


def main():
    realized = set()
    table = Path(__file__).resolve().parent / 'data' / 'brouwer_srg_table_20261118.csv'
    with open(table) as fh:
        for rec in csv.DictReader(fh):
            if rec['status'] in ('+', '!'):
                realized.add(tuple(int(rec[x]) for x in ('v', 'k', 'lam', 'mu')))
    rows = []
    for n, k, lam, mu, r, s, m in smith_sets(RMAX):
        row = (n, k, lam, mu)
        A = thm_main(n, k, lam, mu, r, s, m)
        V = thm_vertex(n, k, lam, mu, r, s, m)
        S2P = second_moment(row)
        subs = subconstituents(n, k, lam, mu, r, s)
        S2S = any(c is not None and second_moment(c) for c in subs)
        u = Fr(k - r, n)
        tight4 = 2 * n == m * (m + 3)
        rows.append((n, k, lam, mu, r, s, str(u), m, feasible(*row),
                     passes_checks(n, k, lam, mu, r, s), A, V, S2P, S2S, tight4, row in realized))
    rows.sort()
    hdr = ('n', 'k', 'lambda', 'mu', 'r', 's', 'u', 'm', 'feasible', 'checks', 'A', 'V', 'S2P',
           'S2S', 'tight4', 'realized')
    if OUT:
        with open(OUT, 'w') as fo:
            fo.write('\t'.join(hdr) + '\n')
            for R in rows:
                fo.write('\t'.join(map(str, R)) + '\n')
    from collections import Counter
    print('Smith sets with r <=', RMAX, ':', len(rows), file=sys.stderr)
    A, V, P, S, C, T4, RL = 10, 11, 12, 13, 9, 14, 15
    print('Theorem 1.1:', Counter(R[A] for R in rows), file=sys.stderr)
    print('Theorem 6.1:', Counter(R[V] for R in rows), file=sys.stderr)
    ex = [R for R in rows if R[A] == 'excl' or R[V] == 'excl']
    ok = [R for R in ex if R[C]]
    print('excluded by 1.1 or 6.1:', len(ex), ' passing the feasibility checks:', len(ok),
          ' of these excluded by Prop 4.4 (graph or subconstituent):', sum(R[P] or R[S] for R in ok),
          ' not:', sum(not (R[P] or R[S]) for R in ok), file=sys.stderr)
    print('Prop 4.4 exclusions missed by 1.1/6.1:',
          sum((R[P] or R[S]) and R[C] and R not in ex for R in rows), file=sys.stderr)
    print('realized:', sum(R[RL] for R in rows), ' realized excluded:',
          sum(R[RL] and (R[A] == 'excl' or R[V] == 'excl' or R[P] or R[S]) for R in rows),
          ' realized failing the checks:', sum(R[RL] and not R[C] for R in rows), file=sys.stderr)
    t4 = [R for R in rows if R[T4]]
    print('tight-4-design (n = m(m+3)/2):', len(t4), ' excluded by 1.1/6.1:',
          sum(R[A] == 'excl' or R[V] == 'excl' for R in t4), file=sys.stderr)


if __name__ == '__main__':
    main()
