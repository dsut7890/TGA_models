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


%morvan's model
%pyrolysis part
T=393:500;
%1.6 k/min 1700 specific heat, 0.1 moisture fraction, heat of pyrolysis
A=0.1; %initial condition

C=1500;
nuChar=8/3;
nuSoot=0;

hpyr=0.418e6;
fudge=10;
beta=fudge*C/(100*hpyr);
%beta=1.6*60*1700*A/(100*hpyr)
%beta=C*(nuChar-nuSoot-1)/(100*hpyr)
gamma=400;
M=A*exp(beta*((T-gamma).^2)/2);


%make the figure
figure(1)
clf
hold on
plot(blackSpruce(:,1),blackSpruce(:,2),'ro')
plot(oak1(:,1),oak1(:,2),'gs')
plot(oak2(:,1),oak2(:,2),'gs')
plot(needles(:,1),needles(:,2),'bd')

%pyrolysis
plot(T,M,'k','linewidth',2)

legend('spruce','oak 1','oak 2', 'pine needles', 'Model - Morvan data','location','southeast')

%moisture loss
plot([393,393],[0,A],'k','linewidth',2)

title('mass loss')
xlabel('T (K)')
ylabel('1-M(T)/M(0)')
set(gca,'fontsize',14)

%char
%figure(2)
%plot(Tchar,mChar,'k','linewidth',2);
