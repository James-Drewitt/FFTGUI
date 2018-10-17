function [G]=FastFT(dat,rho,M,maxq)
%% FastFT(dat,rho,M,Qmax);
%Fourier transform of dat using the FFT algorithm
%      Dr James W E Drewitt
%      Copyright 2018, James W E Drewitt
%      james.drewitt@bristol.ac.uk; james.drewitt@gmail.com
%
% Oringally coded written by James W E Drewitt 01/04/2011, james.drewitt@cnrs-orleans.fr
%
%
% input
%---------
% dat           matrix          column 1 contains Q values. 
%                               column 2 contains S(Q)-1.
% rho           number          number density (A**-3)
% M             number          power of two (N=2^M)          
% Qmax          number          maximum Q-value
% 
% Output
%------------
% results       matrix          column 1 contains r values 
%                               column 2 contains G(r)
%%     
disp('*** Fast Fourier Transform ***');
q=dat(:,1);
sq=dat(:,2);
%
nrows=length(q);
ntrans=(2^M)-1;
nhalf=(ntrans+1)/2;
qstep=q(2)-q(1);
qmax=nhalf*qstep;
nmax=maxq/qstep;
rstep=pi/qmax;
F=zeros(ntrans+1,1);
k=1;
ncut=min(nmax,nrows);
while k<(ncut)
    F(k+1)=sq(k+1)*k*qstep;
    F((ntrans+1)-k+1)=-F(k+1);
    k=k+1;
end
%save 'F.dat' F -ascii
FT=imag(fft(real(F)))/(ntrans);
r=zeros(nhalf,1);
G=zeros(nhalf,2);
for m=1:nhalf
    r(m)=m*rstep;
    G(m+1,1)=r(m);
    G(m,2)=-FT(m)*(qmax/(2*pi*pi*rho*G(m,1)))+1;
end
end