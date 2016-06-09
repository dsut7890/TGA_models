clear
k=1;
rhs=@(t,y) -k*y*exp(-9000/t);

[tout,alpha]=ode45(rhs,[500,900],1);
plot(tout,alpha)
