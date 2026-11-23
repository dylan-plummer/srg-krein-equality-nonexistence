\\ 2-local maximal integral overlattice; check unimodular at 2; characteristic norm mod 4; det odd part mod 4
basisgram(G)={my(KT=qflllgram(G,4),T=KT[2]); T~*G*T};
v2(x)=if(x==0,99,valuation(x,2));
step(B)={my(S=matsnf(B,1),V=S[2],D=S[3],R=#B,gens=List(),t,x,nv);
  for(i=1,R, if(D[i,i]%2==0, listput(gens, V[,i]/2)));
  t=#gens; if(t==0,return(0));
  if(2^t<=2^20,
    print("  step: exhaustive search over ",2^t-1," classes");
    forvec(c=vector(t,i,[0,1]), if(c==vector(t),next);
      x=sum(i=1,t,c[i]*gens[i]); nv=x~*B*x; if(v2(nv)>=0, return(x))),
    print("  step: RANDOM search (2-rank ",t,"); final unimodularity is certified by det");
    for(tr=1,400000, c=vector(t,i,random(2)); if(c==vector(t),next);
      x=sum(i=1,t,c[i]*gens[i]); nv=x~*B*x; if(v2(nv)>=0, return(x))); error("random fail t=",t));
  0};
addvec(B,x)={my(R=#B,den=denominator(x),W); W=mathnf(concat(den*matid(R),den*x))/den; W~*B*W};
run(G,nF)={my(B=basisgram(G),x,d,c,w2,od,M2);
  print("rank ",#B,"  v2(det L) = ",v2(matdet(B)));
  while(1, x=step(B); if(x==0,break); B=addvec(B,x));
  d=matdet(B); print("maximal 2-local integral: v2(det) = ",v2(d), "  (0 = unimodular at 2)");
  \\ characteristic vector: B c == diag(B) mod 2
  M2=Mod(B,2); c=lift(matsolve(M2, Mod(vector(#B,i,B[i,i])~,2)));
  w2=c~*B*c; od=d/2^v2(d);
  print("w^2 mod 4 = ",w2%4,"   n/F mod 4 = ",lift(Mod(numerator(nF),4)/Mod(denominator(nF),4)),
        "   D = odd part of det mod 4 = ",od%4,"   m+1-D mod 4 = ",(#B+1-od)%4)};
