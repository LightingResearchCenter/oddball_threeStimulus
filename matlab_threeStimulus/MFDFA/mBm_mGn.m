function [mBm,mGn]=mBm_mGn(N,Ht)

% Generation of multifractional Brownian motion and multifractional Gaussian noise
%
% [mBm,mGn]=mBm_mGn(N,Ht)
%
% INPUT ARGUMENTS----------------------------------
%
% N:     Sample size
% Ht:    [N,1] vector of the time evolving H(t)
%
% OUTPUT VARIABLES---------------------------------
%
% mBm:   multifractional Brownian motion
% mGn:   multifractional Gaussian noise
%
% EXAMPLE------------------------------------------
%
% N=1024;
% t=1:N;
% Ht=0.5+0.5.*(sin(0.0025.*pi.*t));
% [mBm,mGn]=mBm_mGn(N,Ht);
% figure;
% subplot(311)
% plot(t,Ht);
% ylabel('H');title('Hurst exponent')
% subplot(312)
% plot(t,mBm);
% ylabel('amplitude');title('multifractional Brownian motion')
% subplot(313)
% plot(t,mGn);
% xlabel('time');ylabel('amplitude');title('multifractional Gaussian noise')
%
%--------------------------------------------------
% Written by Espen A. F. Ihlen (espenale@svt.ntnu.no),2009

numb1=10;
numb2=1000;
alpha=2;
N1=numb1*(numb2+N);

mGnSum=zeros(numb1,1);
mGnSumm=zeros(numb1*(numb2-1),1);
mGn=zeros(N,1);

R=randn(N1,1);
for t=1:N,
    for n=1:numb1;
        mGnSum(n)=(n^(Ht(t)-(1/alpha)))*R(1+(numb1*(numb2+t))-n);
    end;
    mGnSum1=sum(mGnSum);
    for nn=1:(numb1*(numb2-1));
        mGnSumm(nn)=(((numb1+nn)^(Ht(t)-(1/alpha)))-nn^(Ht(t)-(1/alpha)))*R(1+(numb1*(numb2-1+t))-nn);
    end;
    mGnSum2=sum(mGnSumm);
    mGn(t)=((numb1^(-Ht(t)))/gamma(Ht(t)-(1/alpha)+1))*(mGnSum1+mGnSum2);
end;

mBm=cumsum(mGn);