clear
%----------------------------------------
%following:
%
%Morvan D. and Dupuy J.L
%"Modeling of Fire Spread Through A Forest
%Fuel Using a Multiphase Formulation"
%Combustion and Flame 127:1981-1994 (2001)
%
%modified to a different reaction model 
%which includes a dependence like A(T)
%Reproduces figure 3 correctly for the 
%two heating rates
%
%----------------------------------------

%all the constants
sigs=700; %surface to volume 
dTdt=10/60; %heating rate

%stoichiometric ratios
nuchar=0.338; %char
nuco2=0.2; %carbon dioxide
nuash=0.033; %ash
nusoot=0.05; %soot
nuo2=8/3; %oxygen

%initial conditions
rho0=270; %initial density
rhog0=1.225; %initial density
alpha0=0.8; %initial occupied volume
mo20=(1-alpha0)*rhog0*0; %initial oxygen
yh2o0=0.08; %initial moisture content
yi0=1-yh2o0; %0.91; %initial dry wood

%temperature vector
T=[300:900];

%constants with error bars
%prefactors
kh2o=6.05e5;
Akpyr0=5.0118e15;
kchar=4.3e2;

%activation energies/R
erh2o=5300;%5556; %h2o FIDDLED!!
erpyr=2.7e4; %pyrolysis
erchar=9000; %char 

%loop over +/- a few percent
perTol=20;
nTol=10;
tol=linspace(-perTol/2,perTol/2,nTol);
tol=tol/100;

%set up the figure
figure(1)
clf
hold on
box on
data01=xlsread('EB_Morvan@10Cpmin.xlsx');
W01=data01(1:100:end,2); % fractional conversion
T01=data01(1:100:end,1); % Temperature
 

%------------------------------------------------
%code begins here
%------------------------------------------------
ics=[alpha0*rho0*yh2o0,alpha0*rho0*yi0,...
    mo20,alpha0,rho0,0,0,0]';

M0=(ics(1)+ics(2)+ics(6));

%the right hand side of the system of equations
%set up for any of the standard ode solvers
%the variables in y
%1: alpha*rho*Y_h2o, 2: alpha*rho*Y_i
%3: oxygen mass, 4: volume fraction
%5: density, 6: alpah*rho*Y_char

A0=0.5*Akpyr0;
A1=1.5*Akpyr0;
T0=300;
T1=900;

Akpyr=@(t) (heaviside(t-T0)-heaviside(t-T1)).*((A1-A0)/(T1-T0).*(t-T0)+A0)+(1-(heaviside(t-T0)-heaviside(t-T1)))*Akpyr0;

%plot([300:900],test([300:900]))
%stop

kpyr=@(a,t) (1/dTdt)*Akpyr(t)*4*(1-(a/M0)).*(-(-(a/M0)-(((a/M0).^2)/2)-(((a/M0).^3)/3)-...
	  (((a/M0).^4)/4)-(((a/M0).^5)/5)-(((a/M0).^6)/6) - ...
	  (((a/M0).^7)/7)-(((a/M0).^8)/8)-(((a/M0).^9)/9) - ...
	  (((a/M0).^10)/10)-(((a/M0).^11)/11)-(((a/M0).^12)/12)- ...
	  (((a/M0).^13)/13)-(((a/M0).^14)/14)-(((a/M0).^15)/15) ...
	  )).^0.75;%(-(log((1-a/M0)))).^(3/4);
%          (((a/M0).^16)/16)-(((a/M0).^17)/17)-(((a/M0).^18)/18) ...
%          (((a/M0).^19)/19)-(((a/M0).^20)/20)-(((a/M0).^21)/21)- ...
%          (((a/M0).^22)/22)-(((a/M0).^23)/23)-(((a/M0).^24)/24)- ...
%          (((a/M0).^25)/25)-(((a/M0).^26)/26)-(((a/M0).^27)/27)- ...
%          (((a/M0).^28)/28)-(((a/M0).^29)/29)-(((a/M0).^30)/30)- ...
%          (((a/M0).^31)/31)-(((a/M0).^32)/32)-(((a/M0).^33)/33) ...



rhs= @(t,y)  [...
	  -(1/dTdt)*((kh2o*y(1)*exp(-erh2o./t)./sqrt(t))).^1;...%moisture mass
	  -(1/dTdt)*kpyr(y(1)+y(2)+y(6),t)*y(2)*exp(-erpyr./t);...%gas mass
	  -(1/dTdt)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs;... %O2 mass
	  -(1/dTdt)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs./(nuo2*y(5));...% alpha
	  (1/dTdt)*((nuchar-nusoot-1)*kpyr(y(1)+y(2)+y(6),t)*y(2)*exp(-erpyr./t) ...
		    -(kh2o*y(1)*exp(-erh2o./t)./sqrt(t)).^1)./y(4);... %rho
	  (1/dTdt)*((nuchar-nusoot)*kpyr(y(1)+y(2)+y(6),t)*y(2)*exp(-erpyr./t) ...
		    -(nuash/nuchar+1)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs/nuo2);...%mass of char
	  (1+nuo2)*(1/dTdt)*kchar*y(3)*exp(-erchar./t)*y(4)*sigs+...
	  nuco2*(1/dTdt)*kpyr(y(1)+y(2)+y(4),t)*y(2)*exp(-erpyr./t);...%co2 fraction
	  (1-nuchar-nuco2)*(1/dTdt)*kpyr(y(1)+y(2)+y(4),t)*y(2)*exp(-erpyr./t);...%co fraction
	  ];

%solve the system, ode23 and ode45 seem stable
[Tnum,sol]=ode23(rhs,T,ics);

%compute the Y values
Yh2o=sol(:,1)./(sol(:,4).*sol(:,5));
Yi=sol(:,2)./(sol(:,4).*sol(:,5));
Yco2=sol(:,6)./(sol(:,4).*sol(:,5));
Yco=sol(:,7)./(sol(:,4).*sol(:,5));
Ychar=1-Yi-Yh2o;
YcharT=sol(:,6)./(sol(:,4).*sol(:,5));

%check conservation
conservationOfMass=max(YcharT-Ychar)

%build the mass loss
mloss=sol(:,1)+sol(:,2)+sol(:,6); %actually mass at time t
mloss=1-mloss/mloss(1);

%-------------------------------------------------
%plot 
%-------------------------------------------------
plot(Tnum,mloss)
plot(T01,W01,'b*')

xlim([300 900]);

title('Mass loss rate','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$1-M/M(0)$','interpreter','latex') 
set(gca,'fontsize',14)

%h=legend('$1-M/M(0)$','Experiment');
%set(h,'position', [0.6738    0.8    0.2137    0.0495]);
%set(h,'interpreter','latex','fontsize',14);
%legend boxoff


