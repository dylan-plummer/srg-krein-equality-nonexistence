"""Audit rational projector-sum coverage of the fifteen BCN files used by the preprint.

Usage in the srg environment: python -B audit_projector_sums.py BCN_DIRECTORY
Printed counts compare the old single-residual enumeration with a split by
irreducible rational factors. Only previously omitted sums need testing here.
"""
import argparse
import itertools
import re
import sys
from collections import Counter
from fractions import Fraction as Q
from math import gcd, lcm
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from distance_arrays import CATALOG, valencies
from frame_test import obstruction


def spectrum(up, down):
    d = len(up)
    older, current = [], [1]
    for i in range(d + 1):
        a = up[0] - (up[i] if i < d else 0) - (down[i-1] if i else 0)
        nxt = current + [0]
        for j, v in enumerate(current):
            nxt[j+1] -= a*v
        if i:
            for j, v in enumerate(older):
                nxt[j+2] -= up[i-1]*down[i-1]*v
        older, current = current, nxt
    p = sp.Poly.from_list(current, gens=sp.Symbol('x'), domain=sp.QQ)
    facts = p.factor_list()[1]
    assert all(e == 1 for _, e in facts)
    return p, [f for f, _ in facts]


def projector(p, f, up, down):
    # CRT polynomial equal to 1 on the roots of f, zero on all other roots.
    other = p.exquo(f)
    indicator = (other * other.invert(f)).rem(p)
    d = len(up)
    v = [Q(0)] * (d+1)
    for coeff in indicator.all_coeffs():
        old = v
        v = []
        for i in range(d+1):
            c = down[i-1] if i else 0
            b = up[i] if i < d else 0
            t = (up[0]-b-c)*old[i] + (c*old[i-1] if i else 0)
            v.append(t + (b*old[i+1] if i<d else 0))
        v[0] += Q(coeff)
    return v


def omitted_sums(p, factors, up, down):
    layers = valencies(up, down)
    n = sum(layers)
    blocks = [projector(p, f, up, down) for f in factors]
    for z in blocks:
        assert sum(h*v*v for h,v in zip(layers,z)) == z[0]
    assert list(map(sum,zip(*blocks))) == [1]+[0]*len(up)
    for choice in itertools.product((0,1), repeat=len(blocks)):
        if len({i for i,f in zip(choice,factors) if f.degree()>1}) < 2:
            continue
        proj = list(map(sum, zip(*(b for b,i in zip(blocks,choice) if i))))
        row_sum = sum(h*v for h,v in zip(layers,proj))
        rank = n*proj[0]
        assert rank.denominator == 1 and row_sum in (0,1)
        assert sum(h*v*v for h,v in zip(layers,proj)) == proj[0]
        scale = lcm(*(v.denominator for v in proj))
        g = [int(scale*v) for v in proj]
        assert gcd(*g) == 1
        cubic = row_sum == 0 and sum(h*v**3 for h,v in zip(layers,g)) == 0
        if row_sum and g[0]%2==0 and scale%2==0:
            yield None
        else:
            yield obstruction(n,int(rank),scale,g,cubic)


def main(directory):
    arrays = set()
    for name in CATALOG:
        for u,d in re.findall(r'i\(([\d,\s]+);([\d,\s]+)\)',(directory/name).read_text()):
            arrays.add((tuple(map(int,u.split(','))),tuple(map(int,d.split(',')))))
    old = full = restricted = tested = 0
    hits = []
    census = Counter()
    for up, down in sorted(arrays):
        p, factors = spectrum(up,down)
        nonlinear = sum(f.degree()>1 for f in factors)
        census[nonlinear] += 1
        rational = len(factors)-nonlinear
        restricted += rational-1
        old += 2**(rational+bool(nonlinear))-2
        full += 2**len(factors)-2
        if nonlinear>1:
            for reason in omitted_sums(p,factors,up,down):
                tested += 1
                if reason:
                    hits.append((up,down,reason))
    assert tested == full-old
    print(f'{len(arrays)} arrays, {restricted} nonprincipal rational eigenvalues')
    print('arrays by number of nonlinear factors',sorted(census.items()))
    print(f'{old} original sums, {full} full sums, {tested} additional, {len(hits)} additional failing sums')
    for hit in hits:
        print(hit)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    main(parser.parse_args().directory)
