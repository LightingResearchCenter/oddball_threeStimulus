function sgates=FTran(signal,N)

% Fourier  transformed surrogates with randomized phases
%
% [sgates]=FTran(signal,N)
%
% INPUT PARAMETERS-----------------------
%
% signal:       Input signal
% N:            Number of surrogates
%
% OUTPUT VARIABLE------------------------
%
% sgates:       A [n,N] matrix of N surrogates with n samples each
%
% EXAMPLE--------------------------------
% 
% N=1024;
% t=1:N;
% Ht=0.5+0.5.*(sin(0.0025*pi.*t));
% [mBm,mGn]=mBm_mGn(N,Ht);
% [sgates]=FTran(mGn,1);
% figure;
% subplot(211);
% plot(t,mGn);
% xlabel('time');ylabel('Amplitude');title('multifractal Gaussian noise');
% subplot(212);
% plot(t,sgates);
% xlabel('time');ylabel('Amplitude');title('surrogate');
%
% ---------------------------------------
% Written by Espen A. F. Ihlen (espen.ihlen@ntnu.no), 2009

signal=signal(:);
M=length(signal);

if nargin<2 | isempty(N)==1
   N=1;
end

Coef=fft(signal);
Amp=abs(Coef); 
Ang=angle(Coef);
i=sqrt(-1);
half=floor(M/2);

for j=1:N
   if rem(M,2)==0
      Angran=rand(half-1,1)*2*pi;
      Ang(2:M)=[Angran' Ang(half+1) -flipud(Angran)'];
      Amp=[Amp(1:half+1);flipud(Amp(2:half))];
   else
      Angran=rand(half,1)*2*pi;
      Ang(2:M)=[Angran -flipud(Angran)];
   end
   sgates(:,j)=Amp.*exp(i*Ang);
   sgates(:,j)=real(ifft(sgates(:,j)));
end

