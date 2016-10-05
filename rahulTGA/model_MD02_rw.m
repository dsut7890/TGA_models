clear
%----------------------------------------
%following:
%
%Morvan D. and Dupuy J.L
%"Modeling of Fire Spread Through A Forest
%Fuel Using a Multiphase Formulation"
%Combustion and Flame 127:1981-1994 (2001)
%
%Reproduces figure 3 correctly for the 
%two heating rates
%
%----------------------------------------

%all the constants
sigs=700; %surface to volume 
dTdt=10/60; %heating rate

%prefactors
kh2o=6.05e5;
Akpyr=5.0118e15;

%stoichiometric ratios
nuchar=0.35; %char
nuco2=0.5; %carbon dioxide
nuash=0.22; %ash
nusoot=0.05; %soot
%nuo2=5/3; %oxygen

%activation energies/R
erh2o=5456; %h2o
erpyr=2.70e4; %pyrolysis
%erchar=10000; %char 

%initial conditions
rho0=270; %initial density
rhog0=1.225; %initial density
alpha0=0.6; %initial occupied volume
yh2o0=0.09; %initial moisture content
yi0=0.91; %initial dry wood
ychar0=0;
%temperature vector
T=[300:900];

%------------------------------------------------
%code begines here
%------------------------------------------------
ics=[alpha0*rho0*yh2o0,alpha0*rho0*yi0,...
    rho0,alpha0*rho0*ychar0,0,0]';

M0=(ics(1)+ics(2)+ics(4));

%the right hand side of the system of equations
%set up for any of the standard ode solvers
%the variables in y
%1: alpha*rho*Y_h2o, 2: alpha*rho*Y_i
%3: oxygen mass, 4: volume fraction
%5: density, 6: alpah*rho*Y_char

kpyr=@(a) (1/dTdt)*Akpyr*4*(1-(a/M0)).*(-(-(a/M0)-(((a/M0).^2)/2)-(((a/M0).^3)/3)-(((a/M0).^4)/4)-(((a/M0).^5)/5)-(((a/M0).^6)/6))).^0.75;%(-(log((1-a/M0)))).^(3/4);

rhs= @(t,y)  [...
          -(1/dTdt)*(((kh2o*y(1)*exp(-erh2o./t))./sqrt(t)));... %moisture mass fraction
          -(1/dTdt)*kpyr(y(1)+y(2)+y(4))*y(2)*exp(-erpyr./t);...%gas ms fraction
          (1/dTdt)*((nuchar-nusoot-1)*kpyr(y(1)+y(2)+y(4))*y(2)*exp(-erpyr./t) ...
	            -((kh2o*y(1)*exp(-erh2o./t))./sqrt(t)).^1);... %oxygen mass
          (1/dTdt)*((nuchar-nusoot)*kpyr(y(1)+y(2)+y(4))*y(2)*exp(-erpyr./t));...%volume fraction
          nuco2*(1/dTdt)*kpyr(y(1)+y(2)+y(4))*y(2)*exp(-erpyr./t);...%co2 fraction
          (1-nuchar-nuco2)*(1/dTdt)*kpyr(y(1)+y(2)+y(4))*y(2)*exp(-erpyr./t);...%co fraction
	  ];

%solve the system, ode23 and ode45 seem stable
[Tnum,sol]=ode45(rhs,T,ics);

%compute the Y values
Yh2o=sol(:,1)./(rho0*alpha0);
Yi=sol(:,2)./(rho0*alpha0);
Yco2=sol(:,5);
Yco=sol(:,6);
Ychar=1-Yi-Yh2o;

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
plot(Tnum,sol(:,4),'b')
title('Quantities','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
legend('alpha*rho*Yh2o','alpha*rho*Yi','alpha*rho*Ychar');
set(gca,'fontsize',14)

figure(3)
clf
hold on
plot(Tnum,sol(:,3),'k')
title('Volume fraction','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$\alpha\, (1/m)$','interpreter','latex') 
set(gca,'fontsize',14)

figure(5)
clf
hold on
%plot(Tnum,Yh2o,'r')
%plot(Tnum,Yi,'g')
plot(Tnum,Ychar,'b')
%plot(Tnum,Yco2,'y')
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

data01=xlsread('EB_Morvan@10Cpmin.xlsx');
xlim([300 900]);
W01=data01(1:50:end,2); % fractional conversion
T01=data01(1:50:end,1); % Temperature
plot(Tnum,mloss,'k',T01,W01,'b*')
saveas(6,'mass loss2','png')

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
plot(ax2,Tnum,mloss,'k',T01,W01,'b*')
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
saveas(7,'output simulatedTGA2','png')

