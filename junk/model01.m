clear

%following combustion and flame 2001
%all the constants
kh2o=6.05e5;
kpyr=3.64e3;
kchar=430;
nuchar=0.338;
nuco2=0.2;
nuash=0.033;
nusoot=0.05;
nuo2=8/3;
erh2o=5956;
erpyr=7250;
echar=9000;
asg=0.5;
sigs=4000; %guessed for pine needles
rhog=1;

%heating rate
dTdt=10/60; 

%get the production of all species
%ompyr=@(as,rhos,ys,Ts) kpyr*as*rhos*ys*exp(-erpyr./Ts);
%omh2o=@(as,rhos,ys,Ts) kh2o*as*rhos*ys*exp(-erh2o./Ts)./sqrt(Ts);
%omchar=@(as,rhos,ys,Ts) kchar*ag*rhog*ys*exp(-erh2o./Ts)*as*sigs./nuo2;

%moisture loss
[th2o,arhoyh2o]=ode45(...
          @(t,y) -(1/dTdt)*((kh2o*y*exp(-erh2o./t)./sqrt(t))).^1,...
          [300,900],0.1);

%loss of the ith species
[ti,arhoyi]=ode45(...
%          @(t,y) -(1/dTdt)*kpyr*y*exp(-erpyr./t),...

          [300,900],0.9);

%mass of O2
%%rhs=@(t,y) [-kchar*y(2)*exp(-erchar./t)*y(1)*sigs,...
%%            -kchar*y(2)*exp(-erchar./t)*y(1)*sigs/rhos

%standard grid
stdT=[300:900];
arhoyh2o=interp1(th2o,arhoyh2o,stdT);
arhoyi=interp1(ti,arhoyi,stdT);


%temporary plot
figure(1)
clf
hold on
plot(stdT,arhoyi,'g')
plot(stdT,arhoyh2o,'b')



stop


%junk?
yh2o=arhoyh2o./arho;
plot(stdT(1:300),yh2o(1:300),'b--')

yh2o=arhoyh2o./arho;
yi=arhoyi./arho;

%make the plot
figure(1)
clf
hold on
plot(T,yh2o)
plot(T,yi)


