function out=kpyr(a,M0,n);

temp=0;
for ii=1:n
   ((a/M0).^ii)/ii
   temp=temp-((a/M0).^ii)/ii;
end

out=-4*(1-(a/M0)).*temp.^(3/4);
