
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

%h=legend('Simulation','Experiment',)

%set(h,'interpreter','latex','fontsize',14);
%set(h,'position', [0.6738    0.5    0.2137    0.0495]);
%set(h,'position', 'northwest');

title('Mass loss rate','interpreter','latex')
xlabel('$T \, (K)$','interpreter','latex') 
ylabel('$1-M/M(0)$','interpreter','latex') 

xlim([300 900]);
legend boxoff

set(gca,'fontsize',14)


