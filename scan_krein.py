#!/usr/bin/env python3
"""Usage: python scan_krein.py NMAX RMAX [OUT.tsv]   (e.g. 10000000 400 hits.tsv)
Independent scan (written from scratch) of ALL srg parameter sets at Krein equality
q^s_ss = 0 (equivalently E(1-2u) = r(r+1)), n <= NMAX, with the characteristic-vector test:

  Odd-norm frame p_i (Gram b(rI-A+uJ), u=a/b), F = bE, f = v2(F).
  (T0) v2(n) >= min(f,4)                         [F w^2 == n mod 16, w^2 in Z_(2)]
  if f <= 2 the maximal integral 2-local lattice is unimodular and w^2 == n/F mod 4:
  (T1) n/F == m (mod 2)                          [D odd]
  (T2) v2(det) even, and n/F == m+1-D (mod 4), when D mod 4 is determined:
       det = b^(m+1) * (primes dividing L) in Q*/Q*^2,  L = gcd(E, numerator(Eu));
       neighbor basis when lambda=0 (det = (br)^(k-1)(br+ak)).
Also records parent/subconstituent feasibility (Cameron-Goethals-Seidel local srgs)."""
import sys, math
from fractions import Fraction as F
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 10**7
RMAX=int(sys.argv[2]) if len(sys.argv)>2 else 400
OUT=sys.argv[3] if len(sys.argv)>3 else None   # output TSV path; omitted = do not write
def v2(x):
    x=F(x)
    if x==0: return 99
    a,b=x.numerator,x.denominator; c=0
    while a%2==0: a//=2; c+=1
    while b%2==0: b//=2; c-=1
    return c
def oddmod4(x):
    x=F(x); a,b=abs(x.numerator),x.denominator
    while a%2==0: a//=2
    while b%2==0: b//=2
    return (a*pow(b,-1,4))%4
def mod4(x):  # x in Z_(2)
    x=F(x); assert v2(x)>=0
    b=x.denominator
    return (x.numerator*pow(b,-1,4))%4
def primes(x):
    out=[];p=2;x=abs(x)
    while p*p<=x:
        if x%p==0:
            out.append(p)
            while x%p==0: x//=p
        p+=1
    if x>1: out.append(x)
    return out
def srg_ok(n,k,l,m):
    if not (0<m<k<n-1 and l>=0) or k*(k-l-1)!=m*(n-k-1): return None
    D=(l-m)**2+4*(k-m); t=math.isqrt(D)
    if t*t!=D: return None
    r=F(l-m+t,2); s=F(l-m-t,2)
    fs=F(k+(n-1)*r, r-s); fr=n-1-fs
    if fs.denominator!=1 or fr.denominator!=1 or fs<=0 or fr<=0: return None
    return r,s,int(fr),int(fs)
def krein_ok(n,k,r,s):
    # both Krein conditions (BvM 1.3.4)
    c1=(r+1)*(k+r+2*r*s)<=(k+r)*(s+1)**2
    c2=(s+1)*(k+s+2*r*s)<=(k+s)*(r+1)**2
    return c1 and c2
def feasible(n,k,l,m,depth=0):
    t=srg_ok(n,k,l,m)
    if t is None: return False
    r,s,fr,fs=t
    if not krein_ok(n,k,r,s): return False
    if n>fr*(fr+3)//2 or n>fs*(fs+3)//2: return False
    return True
def test(n,k,l,m_,r,E):
    s=r-E; u=F(k-r,n); a,b=u.numerator,u.denominator
    m=(k+(n-1)*r)//E
    Fb=b*E; f=v2(Fb)
    if (b*r+a)%2==0:
        # even-norm frame (Theorem 3.1(c)): F w^2 == 0 mod 8; if v2(F)<=1 then m+1-D == 0 mod 4
        if f>=2: return 'n/a(even norm, v2F>=2)',None
        d=E*u; L=math.gcd(E,d.numerator)
        if L%2==1 and ((m+1)*v2(b))%2: return 'EXCL E v2(det) odd',None
        if any(p%4==3 for p in primes(L) if p%2): return 'undetermined D (even norm)',None
        D=pow(oddmod4(b),m+1,4)
        if (m+1-D)%4: return 'EXCL E det (even norm)',(0,(m+1-D)%4)
        return 'pass (even norm)',D
    if v2(n)<min(f,4): return 'EXCL T0',None
    if f>=3: return 'n/a(v2F>=3)',None
    t=mod4(F(n,Fb))
    if (t-m)%2: return 'EXCL T1 parity',t
    d=E*u; L=math.gcd(E,d.numerator)
    bad=[p for p in primes(L) if p%2==1 and p%4==3]
    det_v2_known = (L%2==1)
    if l==0 and m==k:
        disc=F(b*r)**((k-1)%2)*(b*r+a*k)
        if v2(disc)%2: return 'EXCL T2 v2(det) odd (nbr basis)',t
        D=oddmod4(disc)
    else:
        if det_v2_known and ((m+1)*v2(b))%2: return 'EXCL T2 v2(det) odd',t
        if bad: return 'undetermined D (primes %s)'%bad,t
        D=pow(oddmod4(b),m+1,4)
    if (t-(m+1-D))%4: return 'EXCL T2 det',(t,(m+1-D)%4)
    return 'pass',(t,D)
def local_feasible(n,k,l,m,r,E):
    u=F(k-r,n); out=[]
    for order,deg,gap in ((k,l,E*u),(n-k-1,k-m,E*(1-u))):
        if deg==0: out.append('coclique'); continue
        if gap.denominator!=1: return 'local gap nonintegral'
        ls=r-gap; cm=deg+r*ls; cl=cm+r+ls
        if not feasible(order,deg,int(cl),int(cm)) and not (deg==order-1):
            return 'local srg(%d,%d,%s,%s) infeasible'%(order,deg,cl,cm)
        out.append('srg(%d,%d,%s,%s)'%(order,deg,cl,cm))
    return 'ok '+' / '.join(out)
rows=[]
for r in range(1,RMAX+1):
    h=r*(r+1)
    for E in range(h+1, 10**9):
        s=r-E
        P=E-h; Q=2*E          # u=P/Q
        num=Q*(E*r*Q+E*P-h*Q); den=P*(Q-P)
        if den<=0: continue
        # n is increasing in E beyond the minimum; stop when n exceeds NMAX and E large
        nF=F(num,den)
        if nF>NMAX:
            if E>4*h+10: break
            continue
        if nF.denominator!=1: continue
        n=int(nF); kF=r+n*F(P,Q)
        if kF.denominator!=1: continue
        k=int(kF); mu=k+r*s; lam=mu+r+s
        if srg_ok(n,k,lam,mu) is None: continue
        rows.append((n,k,lam,mu,r,E))
rows.sort()
print('Krein-equality (q^s_ss=0) parameter sets with n<=%d, r<=%d: %d'%(NMAX,RMAX,len(rows)))
import collections
C=collections.Counter(); hits=[]
for (n,k,l,m,r,E) in rows:
    v,info=test(n,k,l,m,r,E); C[v.split(' (')[0].split('(')[0]]+=1
    if v.startswith('EXCL'):
        hits.append((n,k,l,m,r,E,v,feasible(n,k,l,m),local_feasible(n,k,l,m,r,E)))
print(C)
if OUT:
  with open(OUT,'w') as fo:
    fo.write('n\tk\tlambda\tmu\tr\tE\ttest\tparent_feasible\tsubconstituents\n')
    for h in hits: fo.write('\t'.join(map(str,h))+'\n')
print('excluded',len(hits),'; parent-feasible & local-ok:',sum(1 for h in hits if h[7] and h[8].startswith('ok')))
for h in hits[:40]: print(*h)
