function [results] = Ipeak(dat)
N=size(dat,1);
r=dat(:,1);
rstep=r(2,1)-r(1,1);
peak=dat(:,2);
intPeak(:,1)=r(:,1);
intPeak(:,2)=0;
for k=2:N
    intPeak(k,2)=intPeak(k-1,2)+dat(k,2).*r(k,1).*r(k,1).*rstep;
end
results=intPeak;