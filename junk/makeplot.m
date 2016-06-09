clear
%load morvan's data
blackSpruce=importdata('black_spruce.csv',',',1);
blackSpruce=blackSpruce.data;
oak1=importdata('oak1.csv',',',1);
oak1=oak1.data;
oak2=importdata('oak2.csv',',',1);
oak2=oak2.data;
needles=importdata('needles.csv',',',1);
needles=needles.data;


%load rahul's data
eucalyptusTimber=importdata('eucalyptus_timber.csv',',',1);
eucalyptusTimber=eucalyptusTimber.data;
eucalyptusBark=importdata('eucalyptus_bark.csv',',',1);
eucalyptusBark=eucalyptusBark.data;
sp=150; %spacing of rahul's data

%morvan's model
%pyrolysis part
T=393:500;
%1.6 k/min 1700 specific heat, 0.1 moisture fraction, heat of pyrolysis
A=0.1; %initial condition

C=1700;
nuChar=8/3;
nuSoot=0;

hpyr=0.418e6;
%beta=C/(100*hpyr);
%beta=1.6*60*1700*A/(100*hpyr)
beta=C*(nuChar-nuSoot-1)/(100*hpyr)
gamma=400;
M=A*exp(beta*((T-gamma).^2)/2);

%rahul pyrolysis
Trahul=393:650;
%1.6 k/min 1700 specific heat, 0.1 moisture fraction, heat of pyrolysis
Arahul=0.08; %initial condition
hpyrRahul=1.8e6; %guessed
beta=1.6*60*1700*Arahul/(100*hpyrRahul);

gamma=400;
Mrahul=Arahul*exp(beta*((Trahul-gamma).^2)/2);


%char part
rhoO2=1;
rhoS=900;
sig=4000;
dTdt=1.6*60;
K1=3*430*sig*rhoO2/(8*dTdt*rhoS);

%numerical solution
%K1=
%[Tchar,mChar]=ode45(@(t,y) y*exp(-9000/t),[500,1100],M(end));

%solve for mass fraction aChar * rhoChar * Ychar
%nuAsh=0.04; 
%nuChar=8/3; %see potire paper
%K2=(nuAsh/nuChar+1);
%[T2,m2]=ode45(@(t,y) K2*K1*y*exp(-9000/t),[500,1100],M(end));


%make the figure
figure(1)
clf
hold on
plot(blackSpruce(:,1),blackSpruce(:,2),'ro')
plot(oak1(:,1),oak1(:,2),'gs')
plot(oak2(:,1),oak2(:,2),'gs')
plot(needles(:,1),needles(:,2),'bd')
plot(eucalyptusBark(1:sp:end,1),eucalyptusBark(1:sp:end,2),'k^')
plot(eucalyptusTimber(1:sp:end,1),eucalyptusTimber(1:sp:end,2),'kv')

%pyrolysis
plot(T,M,'k','linewidth',2)
plot(Trahul,Mrahul,'b','linewidth',2)

legend('spruce','oak 1','oak 2', 'pine needles', 'eucalyptus bark', 'eucalyptus timber', 'Model - Morvan data', 'Model - guessed','location','southeast')

%moisture loss
plot([393,393],[0,A],'k','linewidth',2)

title('mass loss')
xlabel('T (K)')
ylabel('1-M(T)/M(0)')
set(gca,'fontsize',14)

%char
%figure(2)
%plot(Tchar,mChar,'k','linewidth',2);
