function [sgates]=iaaft(signal,N)
% Iterated Amplitude Adjusted Fourier Transformation (IAAFT)
%
% [sgates]=iaaft(signal,N)
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
% [sgates]=iaaft(mGn,1);
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

max=1000;
x=signal(:);
ln=length(x);
amp=abs(fft(x));
sgates=zeros(ln,N);
for n=1:N
    s=randperm(ln);
    sgates(:,n)=x(s);
end

[x,ind]=sort(x);

for n=1:N
    phase_x=angle(fft(sgates(:,n)));
    nn=1;
    conv=0;
    ind_prev=ind;
    while nn<=max && conv==0 
        sgates(:,n)=amp.*exp(phase_x.*1i); 
        sgates(:,n)=real(ifft(sgates(:,n)));
        [sgates(:,n),indx]=sort(sgates(:,n));
        [sgates(:,n),ind_new]=sort(indx);
        sgates(:,n)=x(ind_new);
        if ind_new==ind_prev
            conv=1;
        else
            ind_prev=ind_new;
            nn=nn+1;
        end
        phase_x=angle(fft(sgates(:,n)));
    end    
end
