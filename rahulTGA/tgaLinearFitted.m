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

figure(1)
clf
hold on
box on

%all the constants
sigs=700; %surface to volume 
dTdt=10/60; %heating rate

%load data
data01=xlsread('EB_Morvan@10Cpmin.xlsx');
W01=data01(1:end,2); % fractional conversion
T01=data01(1:end,1); % Temperature

%crop W01 to maximum of 900K 
W01(T01>900)=[];
T01(T01>900)=[];

%pyrolysis range
T0=480;
T1i=650;
DT=5;

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


ctr=0;
n=8;
for ii=1:n
  T0=T0+DT;
  T1=T1i;
  for jj=1:n
    T1=T1+DT;
    %----------------------------
    %get gradient coeffiecent from a fit
    Wt=W01;
    dw=Wt(max(find(T01<=T0)))

    Wt(T01>T1)=[];
    Wt(T01<T0)=[];
    Wt=Wt-dw;
    Tt=T01;
    Tt(T01>T1)=[];
    Tt(T01<T0)=[];

    lna=sum(log(Wt)-2*log(Tt-T0))/length(Tt);
    a=exp(lna);
    a=a*(1+nuchar-nusoot);

    %clf
    %hold on
    %box on
    %plot(T01,W01)
    %plot(Tt,a*(Tt-T0).^2+dw,'r')
    %axis([300,900,0,1])

    %stop

    %------------------------------------------------
    %ode solve code begins here
    %------------------------------------------------
    ics=[yh2o0,yi0,0]';

    M0=(ics(1)+ics(2)+ics(3))*alpha0*rho0;

    %the right hand side of the system of equations
    %set up for any of the standard ode solvers
    %the variables in y
    %1: Y_h2o, 2: Y_i
    %3: Y_char

    %pyr= @(x) (heaviside(x-T0)-heaviside(x-T1)).*(Cpodh.*(1/dTdt)).*(x-T0)/(T1-T0);
    pyr= @(x) 2*(heaviside(x-T0)-heaviside(x-T1)).*a.*(x-T0);

    rhs= @(t,y)  [...
	      -(1/dTdt)*((kh2o*y(1)*exp(-erh2o./t)./sqrt(t))).^1;...%moisture mass
	      -pyr(t);...%mass of substance
	      (nuchar-nusoot).*pyr(t);...%mass of char
	      ];

    %solve the system, ode23 and ode45 seem stable
    [Tnum,sol]=ode23(rhs,T,ics);

    %build the mass loss
    mloss=sol(:,1)+sol(:,2)+sol(:,3); %actually mass at time t
    mloss=1-mloss/mloss(1);

    %compute a correlation coefficient
    Wint=interp1(T01,W01,Tnum);
    Wint=Wint(50:end-50);
    ctr=ctr+1
    err(ctr,:)=[norm(Wint-mloss(50:end-50),2),a,T0,T1];
%-------------------------------------------------
%plot everything for examination
%-------------------------------------------------
    plot(Tnum,mloss,'k',T01,W01,'b*')
    drawnow

  end
end
h=legend('Simulation','Experiment')

set(h,'interpreter','latex','fontsize',14);
%set(h,'position', [0.6738    0.5    0.2137    0.0495]);
set(h,'position', 'northwest');

title('Mass loss rate','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$1-M/M(0)$','interpreter','latex') 

xlim([300 900]);
legend boxoff

set(gca,'fontsize',14)


%find the best fit and plot
[val,ind]=min(err(:,1));
a=err(ind,2);
T0=err(ind,3);
T1=err(ind,4);
%------------------------------------------------
%ode solve code begins here
%------------------------------------------------
ics=[yh2o0,yi0,0]';
M0=(ics(1)+ics(2)+ics(3))*alpha0*rho0;

%the right hand side of the system of equations
%set up for any of the standard ode solvers
%the variables in y
%1: Y_h2o, 2: Y_i
%3: Y_char

%pyr= @(x) (heaviside(x-T0)-heaviside(x-T1)).*(Cpodh.*(1/dTdt)).*(x-T0)/(T1-T0);
pyr= @(x) 2*(heaviside(x-T0)-heaviside(x-T1)).*a.*(x-T0);

rhs= @(t,y)  [...
     -(1/dTdt)*((kh2o*y(1)*exp(-erh2o./t)./sqrt(t))).^1;...%moisture mass
     -pyr(t);...%mass of substance
      (nuchar-nusoot).*pyr(t);...%mass of char
      ];

%solve the system, ode23 and ode45 seem stable
[Tnum,sol]=ode23(rhs,T,ics);

%build the mass loss
mloss=sol(:,1)+sol(:,2)+sol(:,3); %actually mass at time t
mloss=1-mloss/mloss(1);

plot(Tnum,mloss,'r.--'); 
