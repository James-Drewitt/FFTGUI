function [BT]=BFFT(rdat,rho,M)
%% BFFT(dat,rho,M);
%Fourier backtransform transform of real space data using the FFT algorithm
%      Dr James W E Drewitt
%      Copyright 2018, James W E Drewitt
%      james.drewitt@bristol.ac.uk; james.drewitt@gmail.com
%
% Originally written: James W E Drewitt 01/04/2011, james.drewitt@cnrs-orleans.fr
%
%
% input
%---------
% dat           matrix          column 1 contains r values. 
%                               column 2 contains G(r)-1.
% rho           number          number density (A**-3)
% M             number          power of two (N=2^M)          
%  
% Output
%------------
% results       matrix          column 1 contains Q values 
%                               column 2 contains BT(Q)
%%     
disp('*** Back FFT ***');
r=rdat(:,1);
gr=rdat(:,2);
%
nrows=length(r);
ntrans=(2^M)-1;
nhalf=(ntrans+1)/2;
rstep=r(2)-r(1);
rmax=rdat(nrows);
qstep=pi/rmax;
F=zeros(ntrans+1,1);
k=1;
while k<(nrows)
    F(k+1)=gr(k+1)*k*rstep;
    F((ntrans+1)-k+1)=-F(k+1);
    k=k+1;
end
FT=imag(fft(real(F)))/ntrans;
q=zeros(nhalf,1);
BT=zeros(nhalf,2);
for m=1:nhalf
    q(m)=m*qstep;
    BT(m+1,1)=q(m);
    BT(m,2)=-FT(m)*4*pi*rho*rmax/q(m);
end
end