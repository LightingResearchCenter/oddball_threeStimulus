function [mLm,mLn]=mLm_mLn(N,Ht,alpha)

% Generation of multifractional alpha-stable motion and multifractional alpha-stable noise
%
% [mLm,mLn]=mLm_mLn(N,Ht,alpha);
%
% INPUT ARGUMENTS----------------------------------
%
% N:       Sample size
% Ht:      [N,1] vector of the time evolving H(t)
% alpha:   Scaling of the alpha-stable distribution with range [0,2] 
%
% OUTPUT VARIABLES---------------------------------
%
% mLm:     multifractional alpha-stable motion
% mLn:     multifractional alpha-stable noise
%
% EXAMPLE------------------------------------------
%
% N=1024;
% t=1:N;
% Ht=0.5+0.5.*(sin(0.0025.*pi.*t));
% [mLm,mLn]=mLm_mLn(N,Ht,1.7);
% figure;
% subplot(311)
% plot(t,Ht);
% ylabel('H');title('Hurst exponent')
% subplot(312)
% plot(t,mLm);
% ylabel('amplitude');title('multifractional alpha-stable motion')
% subplot(313)
% plot(t,mLn);
% xlabel('time');ylabel('amplitude');title('multifractional alpha-stable noise')
%
%--------------------------------------------------
% Written by Espen A. F. Ihlen (espen.ihlen@ntnu.no),2010

numb1=10;
numb2=1000;
%alpha=2;
N1=numb1*(numb2+N);

mLnSum=zeros(numb1,1);
mLnSumm=zeros(numb1*(numb2-1),1);
mLn=zeros(N,1);
R=alpha_stable_rand(N1,alpha,1);

for t=1:N,
    for n=1:numb1;
        mLnSum(n)=(n^(Ht(t)-(1/alpha)))*R(1+(numb1*(numb2+t))-n);
    end;
    mLnSum1=sum(mLnSum);
    for nn=1:(numb1*(numb2-1));
        mLnSumm(nn)=(((numb1+nn)^(Ht(t)-(1/alpha)))-nn^(Ht(t)-(1/alpha)))*R(1+(numb1*(numb2-1+t))-nn);
    end;
    mLnSum2=sum(mLnSumm);
    mLn(t)=((numb1^(-Ht(t)))/gamma(Ht(t)-(1/alpha)+1))*(mLnSum1+mLnSum2);
end;

mLm=cumsum(mLn);