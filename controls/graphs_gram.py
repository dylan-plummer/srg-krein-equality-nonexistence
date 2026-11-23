#!/usr/bin/env python3
"""Write the scaled centered Gram b(rI-A+uJ) for Higman-Sims and Q^-(5,5) (both checked as srgs).
Used by twoadic.gp for the lattice controls in Section 9 of the paper (see README.md)."""
import itertools, sys
from fractions import Fraction as F
def hs():
    gen=sum(1<<i for i in (11,9,7,6,5,1,0)); blocks=[]
    for msg in range(4096):
        w=0
        for i in range(12):
            if (msg>>i)&1: w^=gen<<i
        w|=(bin(w).count('1')%2)<<23
        if (w&3)==3 and bin(w).count('1')==8: blocks.append(w>>2)
    adj=[set() for _ in range(100)]
    def e(i,j): adj[i].add(j); adj[j].add(i)
    for i in range(22): e(0,1+i)
    for j,B in enumerate(blocks):
        for i in range(22):
            if (B>>i)&1: e(1+i,23+j)
        for l in range(j):
            if not B&blocks[l]: e(23+j,23+l)
    return adj
def qminus(q):
    a,b=next((a,b) for a in range(q) for b in range(1,q) if all((x*x+a*x+b)%q for x in range(q)))
    Q=lambda v:(v[0]*v[1]+v[2]*v[3]+v[4]*v[4]+a*v[4]*v[5]+b*v[5]*v[5])%q
    pts=[v for v in itertools.product(range(q),repeat=6) if any(v) and Q(v)==0 and next(x for x in v if x)==1]
    Bf=lambda u,v:(Q([x+y for x,y in zip(u,v)])-Q(u)-Q(v))%q
    return [set(j for j in range(len(pts)) if j!=i and Bf(pts[i],pts[j])==0) for i in range(len(pts))]
def check(adj):
    n=len(adj); k=len(adj[0]); assert all(len(a)==k for a in adj)
    lm=set()
    for i in range(n):
        for j in range(i+1,n):
            lm.add((j in adj[i], len(adj[i]&adj[j])))
    lam=[c for t,c in lm if t]; mu=[c for t,c in lm if not t]
    assert len(set(lam))==1 and len(set(mu))==1
    return n,k,lam[0],mu[0]
for name,adj,r in (('hs',hs(),2),('q5',qminus(5),4)):
    n,k,l,m=check(adj); u=F(k-r,n); a,b=u.numerator,u.denominator
    G=[[ (b*r+a if i==j else (a-b if j in adj[i] else a)) for j in range(n)] for i in range(n)]
    open(f'G_{name}.gp','w').write('G=['+';'.join(','.join(map(str,row)) for row in G)+'];\n')
    print(name,(n,k,l,m),'u',u,'F',b*(r-(l-m-r)) if False else None)
