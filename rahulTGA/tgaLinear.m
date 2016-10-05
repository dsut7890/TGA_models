clear

%----------------------------------------
%following:
%
%Morvan D. and Dupuy J.L
%Modeling the propagation of a wildfire 
%through a Mediterranean shrub
%using a multiphase formulation
%Combustion and Flame 138:199-210 (2004)
%
%use the linear model for pyrolysis
%to get degradation, keep alpha, sigma, rho 
%constant, and assume Mo2=0
%
%----------------------------------------

%all the constants
sigs=700; %surface to volume 
dTdt=10/60; %heating rate

%pyrolysis range
T0=480;
T1=660;

%load data
data01=xlsread('EB_Morvan@10Cpmin.xlsx');
W01=data01(1:100:end,2); % fractional conversion
T01=data01(1:100:end,1); % Temperature

%estimate from data
Cpodh=3*(1.53e-3-4e-5)*10/((T1-T0))

%prefactors, specific heat, and heat of pyrolysis
kh2o=6.05e5;
%Cp=1700;
%delh=0.58e7;

%stoichiometric ratios
nuchar=0.338; %char
nuco2=0.2; %carbon dioxide
nuash=0.033; %ash
nusoot=0.05; %soot
nuo2=8/3; %oxygen

%activation energies/R
erh2o=5300;%5556; %h2o FIDDLED!!

%initial conditions
rho0=270; %initial density
rhog0=1.225; %initial density
alpha0=0.8; %initial occupied volume
yh2o0=0.08; %initial moisture content
yi0=1-yh2o0; %0.91; %initial dry wood

%temperature vector
T=[300:900];

%------------------------------------------------
%code begines here
%------------------------------------------------
ics=[yh2o0,yi0,0]';

M0=(ics(1)+ics(2)+ics(3))*alpha0*rho0;

%the right hand side of the system of equations
%set up for any of the standard ode solvers
%the variables in y
%1: Y_h2o, 2: Y_i
%3: Y_char

pyr= @(x) (heaviside(x-T0)-heaviside(x-T1)).*(Cpodh.*(1/dTdt)).*(x-T0)/(T1-T0);

rhs= @(t,y)  [...
          -(1/dTdt)*((kh2o*y(1)*exp(-erh2o./t)./sqrt(t))).^1;...%moisture mass
          -(1/dTdt).*pyr(t);...%mass of substance
          (1/dTdt).*(nuchar-nusoot).*pyr(t);...%mass of char
	  ];

%solve the system, ode23 and ode45 seem stable
tic
[Tnum,sol]=ode23(rhs,T,ics);
toc

%build the mass loss
mloss=sol(:,1)+sol(:,2)+sol(:,3); %actually mass at time t
mloss=1-mloss/mloss(1);

%-------------------------------------------------
%plot everything for examination
%-------------------------------------------------
figure(1)
clf
hold on
plot(Tnum,sol(:,1),'r')
plot(Tnum,sol(:,2),'g')
plot(Tnum,sol(:,3),'b')
title('Quantities','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
legend('Yh2o','Yi','Ychar');
set(gca,'fontsize',14)

figure(2)
clf
plot(Tnum,mloss,'k',T01,W01,'b*')
h=legend('Simulation','Experiment')

set(h,'interpreter','latex','fontsize',14);
set(h,'position', [0.6738    0.5    0.2137    0.0495]);

title('Mass loss rate','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$1-M/M(0)$','interpreter','latex') 

xlim([300 900]);
legend boxoff

set(gca,'fontsize',14)
stop

figure(7)
clf
hold on
plot(Tnum,Yh2o,'r')
plot(Tnum,Yi,'g')
plot(Tnum,Ychar,'b')
ylabel('$y$','interpreter','latex') 
set(gca,'fontsize',14)
h=legend('$Y_{h2o}$','$Y_{i}$','$Y_{char}$');
set(h,'interpreter','latex','fontsize',14);
set(h,'position', [0.6738    0.5    0.2137    0.0495]);
legend boxoff

%--------------------------------------------
%build the publication plot
%--------------------------------------------

figure(8)
clf
ax1=gca;
hold on
plot(ax1,Tnum,Yh2o,'r')
plot(ax1,Tnum,Yi,'g')
plot(ax1,Tnum,Ychar,'b')
ylabel('$Y$','interpreter','latex') 
set(gca,'fontsize',14)
ax1.YLim=[0,1];

h=legend('$Y_{h2o}$','$Y_{i}$','$Y_{char}$');
set(h,'interpreter','latex','fontsize',14);
set(h,'position', [0.6738    0.5    0.2137    0.0495]);
legend boxoff

ax2=axes;
plot(ax2,Tnum,mloss,'k',T01,W01,'b*')
ax2.Color='none';
ax2.YLim=[0,1];
ax2.YAxisLocation='right';
linkaxes([ax1,ax2],'x')    
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$1-M/M(0)$','interpreter','latex') 
set(gca,'fontsize',14)
h=legend('$1-M/M(0)$','Experiment');
set(h,'position', [0.6738    0.8    0.2137    0.0495]);
set(h,'interpreter','latex','fontsize',14);
legend boxoff

title('Simulated TGA curves','interpreter','latex')
saveas(7,'output simulatedTGA1','png')

