clear
%----------------------------------------
%following combustion and flame 2004
%the moisture and pyrolysis reactions 
%are replaced by a simple peicewise model
%----------------------------------------

%all the constants
sigs=4000; %surface to volume 
dTdt=10/60; %heating rate

%prefactors
%kh2o=6.05e5;
%kpyr=3.64e3; %replaced by simplified model
kchar=430;

%stoichiometric ratios
nuchar=0.338; %char
nuco2=0.2; %carbon dioxide
nuash=0.033; %ash
nusoot=0.05; %soot
nuo2=8/3; %oxygen

%activation energies
%erh2o=5956; %h2o
%erpyr=7250; %pyrolysis %replaced by a simplified model
erchar=9000; %char 

%latent heats
dhvap=2.25e6;
dhpyr=0.418e6;

%specific heat
C=14000;

%steps
Tv=393;
T0=400;
T1=500;

%initial conditions
rho0=800; %initial density
alpha0=0.5; %initial occupied volume
mo20=0; %initial oxygen
yh2o0=0.1; %initial moisture content
yi0=0.9; %initial dry wood

%temperature vector
T=[300:900];

%------------------------------------------------
%code begines here
%------------------------------------------------
ics=[alpha0*rho0*yh2o0,alpha0*rho0*yi0,...
    mo20,alpha0,rho0,0]';

%the right hand side of the system of equations
%set up for any of the standard ode solvers
%the variables in y
%1: alpha*rho*Y_h2o, 2: alpha*rho*Y_i
%3: oxygen mass, 4: volume fraction
%5: density, 6: alpah*rho*Y_char


%rhs= @(t,y)  [...
%%          -(1/dTdt)*((kh2o*y(1)*exp(-erh2o./t)./sqrt(t))).^1;...
%          -(1/dTdt)*(y(4)*y(5)*C*dTdt)*kroneckerDelta(t-Tv)/dhvap;...
%%          -(1/dTdt)*kpyr*y(2)*exp(-erpyr./t);...
%          -(1/dTdt)*(y(4)*y(5)*C*dTdt/dhpyr)*(t-T0)/(T1-T0)*(heaviside(t-T0)-heaviside(t-T1));...
%          -(1/dTdt)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs;...
%          -(1/dTdt)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs./y(5);...
%          (1/dTdt)*((nuchar-nusoot-1)*(y(4)*y(5)*C*dTdt/dhpyr)*(t-T0)/(T1-T0)*(heaviside(t-T0)-heaviside(t-T1)) ...
%	            -(y(4)*y(5)*C*dTdt)*kroneckerDelta(t-T0)/dhvap)./y(4);...
%          (1/dTdt)*((nuchar-nusoot)*(y(4)*y(5)*C*dTdt/dhpyr)*(t-T0)/(T1-T0)*(heaviside(t-T0)-heaviside(t-T1)) ...
%	            -(nuash/nuchar+1)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs);...
%	  ];

rhs= @(t,y)  [...
          -(1/dTdt)*(y(4)*y(5)*C*dTdt/dhpyr)*(t-T0)/(T1-T0)*(heaviside(t-T0)-heaviside(t-T1));...
          -(1/dTdt)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs;...
          -(1/dTdt)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs./y(5);...
          (1/dTdt)*((nuchar-nusoot-1)*(y(4)*y(5)*C*dTdt/dhpyr)*(t-T0)/(T1-T0)*(heaviside(t-T0)-heaviside(t-T1)) ...
	            )./y(4);...
          (1/dTdt)*((nuchar-nusoot)*(y(4)*y(5)*C*dTdt/dhpyr)*(t-T0)/(T1-T0)*(heaviside(t-T0)-heaviside(t-T1)) ...
	            -(nuash/nuchar+1)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs);...
	  ];



%solve the system, ode23 and ode45 seem stable
[Tnum,sol]=ode45(rhs,T,ics(2:end));

%compute the Y values
Yh2o=ics(1)*heaviside(Tnum-Tv)./(sol(:,3).*sol(:,4));
Yi=sol(:,1)./(sol(:,3).*sol(:,4));
Ychar=sol(:,5)./(sol(:,3).*sol(:,4));

%build the mass loss
%mloss=sol(:,1)+sol(:,2)+sol(:,6);
%mloss=1-mloss/mloss(1);

%-------------------------------------------------
%plot everything for examination
%-------------------------------------------------
%figure(1)
%clf
%hold on
%plot(Tnum,sol(:,2),'g')
%plot(Tnum,sol(:,6),'b')
%title('Quantities','interpreter','latex')
%xlabel('$T \, (K)$','interpreter','latex') 
%legend('alpha*rho*Yh2o','alpha*rho*Yi','alpha*rho*Ychar');
%set(gca,'fontsize',14)

%figure(2)
%clf
%hold on
%plot(Tnum,sol(:,5),'k')
%title('Solid density','interpreter','latex')
%xlabel('$T \, (K)$','interpreter','latex') 
%ylabel('$\rho\, (kg/m^3)$','interpreter','latex') 
%set(gca,'fontsize',14)
%
%figure(3)
%clf
%hold on
%plot(Tnum,sol(:,4),'k')
%title('Volume fraction','interpreter','latex')
%xlabel('$T \, (K)$','interpreter','latex') 
%ylabel('$\alpha\, (1/m)$','interpreter','latex') 
%set(gca,'fontsize',14)
%
%figure(4)
%clf
%hold on
%plot(Tnum,sol(:,3),'k')
%title('Oxygen mass','interpreter','latex')
%xlabel('$T \, (K)$','interpreter','latex') 
%ylabel('$M_{0_2}\, (kg)$','interpreter','latex') 
%set(gca,'fontsize',14)
%
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
%plot(Tnum,mloss,'k')
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
%plot(ax2,Tnum,mloss,'k')
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
saveas(7,'simplifiedmodel','png')
