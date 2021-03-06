clear
%----------------------------------------
%following combustion and flame 2001
%
%
%----------------------------------------

%all the constants
sigs=2000; %surface to volume 
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


erh2o=5956; %h2o
erpyr=8000; %pyrolysis
erchar=9000; %char 

%initial conditions
rho0=800; %initial density
alpha0=1.0; %initial occupied volume
mo20=0; %initial oxygen
yh2o0=0.05; %initial moisture content
yi0=1-yh2o0; %initial dry wood

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
rhs= @(t,y)  [...
          -(1/dTdt)*((kh2o*y(1)*exp(-erh2o./t)./sqrt(t))).^1;...
          -(1/dTdt)*kpyr*y(2)*exp(-erpyr./t);...
          -(1/dTdt)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs;...
          -(1/dTdt)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs./y(5);...
          (1/dTdt)*((nuchar-nusoot-1)*kpyr*y(2)*exp(-erpyr./t) ...
	            -(kh2o*y(1)*exp(-erh2o./t)./sqrt(t)).^1)./y(4);...
          (1/dTdt)*((nuchar-nusoot)*kpyr*y(2)*exp(-erpyr./t) ...
	            -(nuash/nuchar+1)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs);...
	  ];

%solve the system, ode23 and ode45 seem stable
[Tnum,sol]=ode45(rhs,T,ics);

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

%load data
ax2=axes;
hold on
plot(ax2,Tnum,mloss,'k')
ax2.Color='none';
ax2.YLim=[0,1];
ax2.YAxisLocation='right';
linkaxes([ax1,ax2],'x')    
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$1-M/M(0)$','interpreter','latex') 
set(gca,'fontsize',14)

%put rahul's data on the plot
eucalyptusTimber=importdata('timbertest.csv',',',1);
eucalyptusTimber=eucalyptusTimber.data;
sp=150; %spacing of rahul's data
plot(eucalyptusTimber(1:sp:end,1),eucalyptusTimber(1:sp:end,2),'kv')

h=legend('$1-M/M(0)$','Eucalyptus');
set(h,'position', [0.6738    0.8    0.2137    0.0495]);
set(h,'interpreter','latex','fontsize',14);
legend boxoff

title('Simulated and experimental TGA curves: Eucalyptus timber','interpreter','latex')
saveas(1,'overlayed','png')


