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



%data
T=400:500;
%1.6 K/min 1700 specific heat, 0.1 fraction, heat of pyrolysis
beta=1.6*60*1700*(0.1)/(100*0.418e6);

gamma=400;
A=0.1;
M=A*exp(beta*((T-gamma).^2)/2);

figure(1)
clf
hold on
plot(T,M)
