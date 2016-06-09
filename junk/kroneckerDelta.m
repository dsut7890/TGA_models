function f=kroneckerDelta(m)

ind=find(m==0);
f=0*m;
f(ind)=1;
