"""Checks for Section 5.3 (second subconstituent of Kre(t)); srg conda env, about 1 minute.

1. Two expansions of Phi(G) for triangle-free srgs (cell partition, as in the paper, and
   homomorphism counts) agree with brute force on
   Petersen, Clebsch, Gewirtz, M22 (77) and Higman-Sims.
2. The identities of the proof of Theorem 5.1 hold in Q(t).
3. Control t=2: from a shuffled srg(77,16,0,4) alone, the Norton algebra recovers 22 points
   and a 3-(22,6,1) design whose disjointness graph is the input graph.
"""

import itertools, numpy as np, sympy as sp
from fractions import Fraction as Fr

PAIRS=list(itertools.combinations(range(4),2))
def phi_formula(n,k,mu,a,b):
    """Phi(G), G=(1-a)I+(b-a)A+aJ, by expanding each of the 6 entries into I/A/J parts
    and counting homomorphisms into a triangle-free srg(n,k,0,mu)."""
    y=b-a; ell=n-1-k; total=0
    for lab in itertools.product('IAJ',repeat=6):
        # union-find on I-edges
        par=list(range(4))
        def f(x):
            while par[x]!=x: x=par[x]
            return x
        for (p,q),L in zip(PAIRS,lab):
            if L=='I': par[f(p)]=f(q)
        cls=sorted({f(x) for x in range(4)})
        edges=set(); bad=False
        for (p,q),L in zip(PAIRS,lab):
            if L=='A':
                P,Q=f(p),f(q)
                if P==Q: bad=True;break
                edges.add(frozenset((P,Q)))
        if bad: continue
        V=len(cls); E=len(edges)
        # triangle check
        tri=any(frozenset((x,y_)) in edges and frozenset((y_,z)) in edges and frozenset((x,z)) in edges
                for x,y_,z in itertools.combinations(cls,3))
        if tri: continue
        # components
        comp={c:c for c in cls}
        def g(x):
            while comp[x]!=x: x=comp[x]
            return x
        for e in edges:
            x,z=tuple(e); comp[g(x)]=g(z)
        C=len({g(c) for c in cls})
        if E==V-C: hom=n**C*k**E           # forest
        else:
            assert V==4 and E==4            # must be C4
            hom=n*(k*k+ell*mu*mu)
        w=(1-a)**lab.count('I')*y**lab.count('A')*a**lab.count('J')
        total+=w*hom
    return total

def brute(Adj,a,b):
    n=len(Adj); G=(1-a)*np.eye(n)+(b-a)*Adj+a*np.ones((n,n))
    W=G[:,None,:]*G[None,:,:]
    return float(np.sum(G*np.einsum('ijk,kl,ijl->ij',W,G,W)))

# --- graphs
def golay_octads():
    g=[1,0,1,0,1,1,1,0,0,0,1,1]  # x^11+x^10+x^6+x^5+x^4+x^2+1, low->high
    gens=[]
    for s in range(12):
        v=[0]*23
        for i,c in enumerate(g):
            if c: v[i+s]=1
        v.append(sum(v)%2); gens.append(int(''.join(map(str,v)),2))
    code={0}
    for gg in gens: code|={c^gg for c in code}
    assert len(code)==4096
    return [c for c in code if bin(c).count('1')==8]
oct_=golay_octads(); assert len(oct_)==759
pts=lambda c:{i for i in range(24) if c>>i&1}
blocks=[pts(c)-{22,23} for c in oct_ if {22,23}<=pts(c)]
assert len(blocks)==77 and all(len(B)==6 for B in blocks)
def disj(Bs): return np.array([[1.0 if (i!=j and not (Bs[i]&Bs[j])) else 0.0 for j in range(len(Bs))] for i in range(len(Bs))])
A77=disj(blocks)
A56=disj([B for B in blocks if 0 not in B])
# HS from the 3-(22,6,1) design
n=100; AHS=np.zeros((n,n))
for p in range(22): AHS[0,1+p]=AHS[1+p,0]=1
for i,B in enumerate(blocks):
    for p in B: AHS[1+p,23+i]=AHS[23+i,1+p]=1
AHS[23:,23:]=A77
pet=list(itertools.combinations(range(5),2))
APet=np.array([[1.0 if not set(x)&set(y) else 0.0 for y in pet] for x in pet])
cv=list(itertools.product([0,1],repeat=4))
ACl=np.array([[1.0 if sum(p!=q for p,q in zip(x,y)) in (1,4) else 0.0 for y in cv] for x in cv])
def params(A):
    n=len(A); k=int(A[0].sum()); A2=A@A
    lam={int(A2[i,j]) for i in range(n) for j in range(n) if A[i,j]}
    mu={int(A2[i,j]) for i in range(n) for j in range(n) if i!=j and not A[i,j]}
    return n,k,lam,mu
for name,A in [('Petersen',APet),('Clebsch',ACl),('Gewirtz',A56),('M22-77',A77),('HS',AHS)]:
    n,k,lam,mu=params(A); assert lam=={0} and len(mu)==1; mu=mu.pop()
    for a,b in [(0.13,-0.31),(-0.2,0.45)]:
        F=phi_formula(n,k,mu,a,b); B=brute(A,a,b)
        assert abs(F-B)<1e-7*max(1,abs(B)),(name,F,B)
    print('formula = brute force on',name,(n,k,0,mu))


# --- cell-partition evaluation of Phi (Theorem 5.1, Step 1)
def phi_cells(n,k,mu,w0,w1):
    kb=n-1-k; y=w1-w0
    def Qval(cells,E):
        S1=sum(s*w for s,w in cells); S2=sum(s*w*w for s,w in cells)
        SA=sum(E[(X,Y)]*cells[X][1]*cells[Y][1] for (X,Y) in E)
        return (1-w0)*S2+y*SA+w0*S1**2
    def complete(cells,E):
        # fill rest-rest from degree count: total ordered pairs from R = size_R*k
        R=len(cells)-1; out=sum(E.get((R,Y),0) for Y in range(len(cells)) if Y!=R)
        E[(R,R)]=cells[R][0]*k-out
        return E
    def sym(E):
        F=dict(E)
        for (X,Y),v in E.items(): F[(Y,X)]=v
        return F
    # i=j: cells {i}, N=Gamma(i), R
    c=[(1,1),(k,w1**2),(kb,w0**2)]
    E=sym({(0,1):k,(1,2):k*(k-1)}); E[(1,1)]=0; E=complete(c,E); Qeq=Qval(c,E)
    # i~j: {i},{j},Ai=Gamma(i)-{j},Aj,R
    c=[(1,w1),(1,w1),(k-1,w1*w0),(k-1,w1*w0),(n-2*k,w0**2)]
    E=sym({(0,1):1,(0,2):k-1,(1,3):k-1,(2,3):(k-1)*(mu-1),(2,4):(k-1)*(k-mu),(3,4):(k-1)*(k-mu)})
    E[(2,2)]=0;E[(3,3)]=0; E=complete(c,E); Qadj=Qval(c,E)
    # i!~j: {i},{j},C,Ai,Aj,R
    c=[(1,w0),(1,w0),(mu,w1**2),(k-mu,w1*w0),(k-mu,w1*w0),(n-2-2*k+mu,w0**2)]
    E=sym({(0,2):mu,(0,3):k-mu,(1,2):mu,(1,4):k-mu,(2,5):mu*(k-2),(3,4):(k-mu)*mu,(3,5):(k-mu)*(k-1-mu),(4,5):(k-mu)*(k-1-mu)})
    for X in (2,3,4): E[(X,X)]=0
    E=complete(c,E); Qnon=Qval(c,E)
    return n*(Qeq+k*w1*Qadj+kb*w0*Qnon), (Qeq,Qadj,Qnon)
for name,A in [('Petersen',APet),('Clebsch',ACl),('Gewirtz',A56),('M22-77',A77),('HS',AHS)]:
    n_,k_,lam_,mu_=params(A); mu_=mu_.pop()
    for a_,b_ in [(0.13,-0.31),(-0.2,0.45)]:
        B_=brute(A,a_,b_); assert abs(float(phi_cells(n_,k_,mu_,a_,b_)[0])-B_)<1e-7*max(1,abs(B_))
    print('cell partition = brute force on',name)

# --- identities in Q(t)
t=sp.symbols('t',positive=True)
v=t*(t**2+3*t+1); K=t*(t+1); n=(t**2+2*t-1)*(t**2+3*t+1); k=t**2*(t+2); mu=t**2; ell=n-1-k
r,s=t,-K
assert sp.expand((sp.Symbol('x')-r)*(sp.Symbol('x')-s)-(sp.Symbol('x')**2+mu*sp.Symbol('x')-(k-mu)))==0
g=sp.factor((k+(n-1)*r)/(r-s)); assert sp.expand(g-(v-1))==0
w1=s/k; w0=-(1+s)/ell
assert sp.simplify(w1+(t+1)/(t*(t+2)))==0 and sp.simplify(w0-1/((t+1)*(t+2)))==0
M=1+k*w1**3+ell*w0**3
assert sp.simplify(M-(t**2+t-1)*(t**2+3*t+1)/(t*(t+1)**2*(t+2)**2))==0
assert sp.simplify(phi_formula(n,k,mu,w0,w1)-n**2*(g-2)*M**2/(g*(g-1)))==0
Pc,(Qe,Qa,Qn)=phi_cells(n,k,mu,w0,w1)
assert sp.simplify(Pc-n**2*(g-2)*M**2/(g*(g-1)))==0 and sp.simplify(Qe-M**2)==0
assert sp.simplify(Pc-(t**2+t-1)*(t**2+2*t-1)*(t**2+3*t+1)**4*(t**3+3*t**2+t-3)/(t**2*(t+1)**5*(t+2)**5))==0
assert sp.simplify(Qa-2*(t**2+t-1)*(t**2+3*t+1)**2/(t**3*(t+1)**2*(t+2)**5))==0
assert sp.simplify(Qn-(t**2+3*t+1)**2*(t**5+7*t**4+15*t**3+11*t**2+t+1)/(t**2*(t+1)**5*(t+2)**5))==0
tau=(g/n)**2*M; kap=tau/(g-1); rho=tau/(v-2*K)
assert sp.simplify(rho**2*K*(v-K)-kap*g/n)==0
for w,size in ((w1,0),(w0,t)): assert sp.simplify((K**2+w*K*(v-K))/v-size)==0
assert sp.simplify(n-v*(v-1)/K)==0 and sp.simplify(n*K*(K-1)*(K-2)/(v*(v-1)*(v-2))-(t-1))==0
print('identities of Theorem 5.1 hold in Q(t); tau =',sp.factor(tau))

# --- control t=2
A=A77.copy(); rng=np.random.default_rng(1); perm=rng.permutation(77); A=A[perm][:,perm]
w,V=np.linalg.eigh(A); U=V[:,np.isclose(w,-6)]; f=U.shape[1]; print('f =',f)
P=U@U.T; n=77
tau=np.sum(P**3)/f
W=P[:,None,:]*P[None,:,:]; Phi=np.sum(P*np.einsum('ijk,kl,ijl->ij',W,P,W))
print('Phi =',Phi,' f(f-2)tau^2/(f-1) =',f*(f-2)*tau**2/(f-1))
# tensor T in orthonormal basis U: T_abcd = <e_a*e_b, e_c*e_d>, e_a*e_b = P(U[:,a]∘U[:,b])
Y=np.einsum('ia,ib->abi',U,U)           # e_a∘e_b in R^n
Z=np.einsum('abi,ic->abc',Y,U)          # coords of P(e_a∘e_b) in basis U
T=np.einsum('abx,cdx->abcd',Z,Z)
D=T-T.transpose(0,2,1,3); I=np.eye(f)
S=np.einsum('ac,bd->abcd',I,I)-np.einsum('ab,cd->abcd',I,I)
K=tau/(f-1); print('max|D-KS| =',np.abs(D-K*S).max())
# algebra: on U+Re, x.y = x*y + K<x,y>e; find idempotents by diagonalizing a generic multiplication operator
c=rng.standard_normal(f); L=np.zeros((f+1,f+1))
L[:f,:f]=np.einsum('c,abc->ab',c,Z)   # x*e_b component: (c*e_b)_a
L[f,:f]=K*c                         # e-component of c.e_b is K<c,e_b>
L[:f,f]=c                           # c.e = c
ev,EV=np.linalg.eig(L); print('distinct eigenvalues:',len(set(np.round(ev.real,8))))
# coordinates of u_i (= U^T e_i) in eigenbasis -> two-valued?
EVi=np.linalg.inv(EV)
def Lmat(c):
    L=np.zeros((f+1,f+1)); L[:f,:f]=np.einsum('c,abc->ab',c,Z); L[f,:f]=K*c; L[:f,f]=c; return L
X=np.array([np.diag(EVi@Lmat(U[i])@EV).real for i in range(77)]).T
blocks=[]
for i in range(77):
    col=X[:,i]; vals=np.unique(np.round(col,6)); assert len(vals)==2
    blocks.append(frozenset(np.where(col>0)[0]))
print('block sizes',{len(B) for B in blocks},' distinct',len(set(blocks)))
cnt={}
for S3 in itertools.combinations(range(22),3):
    cnt[sum(1 for B in blocks if set(S3)<=B)]=1
print('3-subset multiplicities',set(cnt))
ok=all((A[i,j]==1)==(len(blocks[i]&blocks[j])==0) for i in range(77) for j in range(77) if i!=j)
print('adjacency = disjointness:',ok)
