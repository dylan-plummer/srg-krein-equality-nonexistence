"""List the failing rational projector sums (Section 8) for the nine further BCN arrays: for each array, every
nonempty proper set of irreducible factors of the characteristic polynomial (principal idempotent excluded), the sum
of the corresponding primitive idempotents, and the verdict of Theorem 3.1 + Lemma 3.2."""
import itertools, sys
from fractions import Fraction as Fr
from math import gcd
import sympy as sp
sys.path.insert(0,'.')
from frame_check import criterion
ARRAYS=[((39,30,4),(1,5,36)),((25,24,2,1),(1,2,24,25)),((70,69,56,10),(1,14,60,70)),((115,96,20,1),(1,20,96,115)),
        ((15,14,14,11),(1,1,4,15)),((594,550,60,1),(1,30,550,594)),((285,280,6,1),(1,6,280,285)),
        ((165,160,10,1),(1,2,160,165)),((169,168,2,1),(1,2,168,169))]
x=sp.symbols('x')
for b,c in ARRAYS:
    d=len(b); B=list(b)+[0]; C=[0]+list(c); k0=b[0]
    k=[1]
    for i in range(d): k.append(k[-1]*b[i]//c[i])
    n=sum(k)
    T=sp.zeros(d+1,d+1)
    for i in range(d+1):
        T[i,i]=k0-B[i]-C[i]
        if i<d: T[i,i+1]=B[i]
        if i>0: T[i,i-1]=C[i]
    facs=[f for f,_ in sp.factor_list(T.charpoly(x).as_expr())[1] if sp.Poly(f,x).eval(k0)!=0]
    def col(theta):
        u=[sp.Integer(1),theta/k0]
        for i in range(1,d):
            a=k0-B[i]-C[i]; u.append(sp.expand(((theta-a)*u[i]-C[i]*u[i-1])/B[i]))
        m=sp.simplify(n/sum(k[i]*u[i]**2 for i in range(d+1)))
        return [m*ui for ui in u]
    fails=[]
    for r in range(1,len(facs)+1):
        for S in itertools.combinations(range(len(facs)),r):
            if r==len(facs): continue
            tot=[sp.Integer(0)]*(d+1); labels=[]
            for j in S:
                rts=sp.Poly(facs[j],x).all_roots(); labels.append(str(facs[j]) if len(rts)>1 else str(rts[0]))
                for t in rts:
                    cl=col(t); tot=[a+bb for a,bb in zip(tot,cl)]
            tot=[sp.nsimplify(sp.simplify(z)) for z in tot]
            fr=[Fr(int(sp.numer(z)),int(sp.denom(z))) for z in tot]
            L=1
            for q in fr: L=L*q.denominator//gcd(L,q.denominator)
            g=[int(q*L) for q in fr]; G=0
            for z in g: G=gcd(G,z)
            g=[z//G for z in g]
            if g[0]<0: g=[-z for z in g]
            F=Fr(sum(k[i]*g[i]**2 for i in range(d+1)),g[0]); m=Fr(n*g[0])/F
            if F.denominator!=1 or m.denominator!=1: continue
            cub=sum(k[i]*g[i]**3 for i in range(d+1))==0
            v=criterion(n,int(m),int(F),g,cub)
            if v: fails.append(('+'.join(labels),g,int(F),int(m),cub,v))
    print(f'{n}\t{{{",".join(map(str,b))};{",".join(map(str,c))}}}')
    for f in fails: print('   E sum over eigenvalue(s)',f[0],' gamma',f[1],' F',f[2],' m',f[3],' cubic',f[4],' ->',f[5])
