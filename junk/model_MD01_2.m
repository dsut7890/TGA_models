clear
%----------------------------------------
%following combustion and flame 2001
%
%
%----------------------------------------

%all the constants
sigs=4000; %surface to volume 
dTdt=10/60; %heating rate

%prefactors
kh2o=6.05e5;
kpyr=3.64e3;
kchar=430;

%stoichiometric ratios
nuchar=0.338; %char
nuco2=0.2; %carbon dioxide
nuash=0.033; %ash
nusoot=0.05; %soot
nuo2=8/3; %oxygen

%activation energies
erh2o=5956; %h2o
erpyr=7250; %pyrolysis
erchar=9000; %char 

%initial conditions
rho0=800; %initial density
alpha0=0.5; %initial occupied volume
mo20=0; %initial oxygen
mco0=0; %carbon monoxide
mco20=0;%carbon dioxide
mh2o0=0;%water vapour
yh2o0=0.1; %initial moisture content
yi0=0.9; %initial dry wood

%temperature vector
T=[300:900];

%------------------------------------------------
%code begines here
%------------------------------------------------
ics=[alpha0*rho0*yh2o0,alpha0*rho0*yi0,...
    mo20,alpha0,rho0,0]';%,mco0,mco20,mh2o0]';

M0=ics(1)+ics(2)+mo20+mco0+mco20+mh2o0;

%the right hand side of the system of equations
%set up for any of the standard ode solvers
%the variables in y
%1: alpha*rho*Y_h2o, 2: alpha*rho*Y_i
%3: oxygen mass, 4: volume fraction
%5: density, 6: alpah*rho*Y_char

%7: mass co, 8: mass co2
%9: mass h2o


rhs= @(t,y)  [...
          -(1/dTdt)*((kh2o*y(1)*exp(-erh2o./t)./sqrt(t))).^1;...
          -(1/dTdt)*kpyr*y(2)*exp(-erpyr./t);...
          -(1/dTdt)*kchar*( M0-y(1)-y(2)-y(6) )*exp(-erchar./t)*y(4)*sigs;...
          -(1/dTdt)*kchar*( M0-y(1)-y(2)-y(6)  )*exp(-erchar./t)*y(4)*sigs./y(5);...
          (1/dTdt)*((nuchar-nusoot-1)*kpyr*y(2)*exp(-erpyr./t) ...
	            -(kh2o*y(1)*exp(-erh2o./t)./sqrt(t)).^1)./y(4);...
          (1/dTdt)*((nuchar-nusoot)*kpyr*y(2)*exp(-erpyr./t) ...
	            -(nuash/nuchar+1)*kchar*(  M0-y(1)-y(2)-y(6) )*exp(-erchar./t)*y(4)*sigs);...
%	  (1/dTdt)*( (1-nuchar-nuco2)*kpyr*y(2)*exp(-erpyr./t));...
%          (1/dTdt)*((1+nuo2)*kchar*( M0-y(1)+y(2)+y(6) )*exp(-erchar./t)*y(4)*sigs ...
%	            +nuco2*kpyr*y(2)*exp(-erpyr./t));...
%          (1/dTdt)*((kh2o*y(1)*exp(-erh2o./t)./sqrt(t))).^1;...
	  ];

%solve the system, ode23 and ode45 seem stable
[Tnum,sol]=ode15s(rhs,T,ics);

%compute the Y values
Yh2o=sol(:,1)./(sol(:,4).*sol(:,5));
Yi=sol(:,2)./(sol(:,4).*sol(:,5));
Ychar=sol(:,6)./(sol(:,4).*sol(:,5));

%build the mass loss
mloss=sol(:,1)+sol(:,2)+sol(:,6);
mloss=1-mloss/mloss(1);

%-------------------------------------------------
%plot everything for examination
%-------------------------------------------------
figure(1)
clf
hold on
plot(Tnum,sol(:,1),'r')
plot(Tnum,sol(:,2),'g')
plot(Tnum,sol(:,6),'b')
title('Quantities','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
legend('alpha*rho*Yh2o','alpha*rho*Yi','alpha*rho*Ychar');
set(gca,'fontsize',14)

figure(2)
clf
hold on
plot(Tnum,sol(:,5),'k')
title('Solid density','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$\rho\, (kg/m^3)$','interpreter','latex') 
set(gca,'fontsize',14)

figure(3)
clf
hold on
plot(Tnum,sol(:,4),'k')
title('Volume fraction','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$\alpha\, (1/m)$','interpreter','latex') 
set(gca,'fontsize',14)

figure(4)
clf
hold on
plot(Tnum,sol(:,3),'k')
title('Oxygen mass','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$M_{0_2}\, (kg)$','interpreter','latex') 
set(gca,'fontsize',14)

figure(5)
clf
hold on
plot(Tnum,Yh2o,'r')
plot(Tnum,Yi,'g')
plot(Tnum,Ychar,'b')
title('Mass fraction','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$Y$','interpreter','latex') 
set(gca,'fontsize',14)

figure(6)
clf
hold on
plot(Tnum,mloss,'k')
title('Mass loss rate','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$1-M/M(0)$','interpreter','latex') 
set(gca,'fontsize',14)

%--------------------------------------------
%build the publication plot
%--------------------------------------------

figure(7)
clf
ax1=gca;
hold on
plot(ax1,Tnum,Yh2o,'r')
plot(ax1,Tnum,Yi,'g')
plot(ax1,Tnum,Ychar,'b')
ylabel('$Y$','interpreter','latex') 
set(gca,'fontsize',14)
ax1.YLim=[0,1];

h=legend('$Y_{H2O}$','$Y_{i}$','$Y_{Char}$');
set(h,'interpreter','latex','fontsize',14);
set(h,'position', [0.6738    0.5    0.2137    0.0495]);
legend boxoff

ax2=axes;
plot(ax2,Tnum,mloss,'k')
ax2.Color='none';
ax2.YLim=[0,1];
ax2.YAxisLocation='right';
linkaxes([ax1,ax2],'x')    
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$1-M/M(0)$','interpreter','latex') 
set(gca,'fontsize',14)
h=legend('$1-M/M(0)$');
set(h,'position', [0.6738    0.8    0.2137    0.0495]);
set(h,'interpreter','latex','fontsize',14);
legend boxoff

title('Simulated TGA curves','interpreter','latex')
saveas(7,'simulatedTGA','png')
