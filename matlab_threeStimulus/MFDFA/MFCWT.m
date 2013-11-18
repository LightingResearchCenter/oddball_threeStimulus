function [s,q,Hq,h,Dh,logFq]=MFCWT(signal,scmin,scmax,ressc,qmin,qmax,qres)

% Multifractal continuous wavelet transformation (MFCWT)
%
% function [s,q,Hq,h,Dh,logFq]=MFCWT(signal,scmin,scmax,ressc,qmin,qmax,qres)
%
% INPUT PARAMETERS---------------------------------------------
%
% signal:       input signal
% scmin:        lower bound of the window size s
% scmax:        upper bound of the window size s
% ressc:        number of elements in s 
% qmin          lower bound of q
% qmax          upper bound of q
% qres          number of elements in q
%
% OUTPUT VARIABLES--------------------------------------------
%
% s:            scale
% q:            q-orders
% Hq:           q-generalized Hurst exponent
% h:            Hölder exponent
% Dh:           Multifractal spectrum
% logFq:        Logarithm of the scaling function F(s,q)
%
% EXAMPLE-----------------------------------------------------
%
% T=40.96 ; rmin=0.02 ; Dt=rmin/2 ; Law=0 ; ParamLaw1=0.2; ParamLaw2=0; C=1 ; q1=(-5:5)/2;
% [Qr1,time,TtheoQ,q] = IDCnoise_epl(T,rmin,Dt,Law,ParamLaw1,C,q1); 
% [Qr2,time,TtheoQ,q] = IDCnoise_epl(T,rmin,Dt,Law,ParamLaw2,C,q1); 
% H=0.75 ; oversamprate = 16 ; printout = 0;
% [Ar,VH1,BH] = synthArVH(H,Qr1,Dt,oversamprate,printout) ;
% [Ar,VH2,BH] = synthArVH(H,Qr2,Dt,oversamprate,printout) ;
% [s,q,Hq1,h1,Dh1,logFq1]=MFCWT(diff(VH1),10,200,40,0.1,3,30);
% [s,q,Hq2,h2,Dh2,logFq2]=MFCWT(diff(VH2),10,200,40,0.1,3,30);
% figure;
% subplot(311)
% plot(1:4096,diff(VH1),'b-',1:4096,diff(VH2),'r-');xlabel('time');
% ylabel('amplitude'); title('multiplicative cascading noise \DeltaB_H(A(t))and fractional Gaussian noise \DeltaB_H(t)');
% legend('\DeltaB_H(A(t))','\DeltaB_H(t)');
% subplot(323);
% plot(log2(s),logFq1(:,1:5:end),'b-');hold all
% plot(log2(s),logFq2(:,1:5:end),'r-');
% xlabel('log_2(\Deltat)');ylabel('log_2(F_\Delta_t(q))');title('log-scaling function')
% subplot(324)
% plot(q,Hq1,'bo-',q,Hq2,'ro-');xlabel('q'); ylabel('H(q)'); title('q-generalized Hurst exponent');
% legend('\DeltaB_H(A(t))','\DeltaB_H(t)');
% subplot(325)
% plot(h1, Dh1,'bo-',h2, Dh2,'ro-');xlabel('h(q)');ylabel('D(h)');title('multifractal half-spectrum');
% legend('\DeltaB_H(A(t))','\DeltaB_H(t)');
%
% ---------------------------------------------------------------
% Written by Espen A. F. Ihlen (espen.ihlen@ntnu.no), 2009

scale=2.^linspace(log2(scmin/6.5),log2(scmax/6.5),ressc);
q=linspace(qmin,qmax,qres);
signal=signal-mean(signal);
wt = cwt(signal,scale,'morl');
s=2.^linspace(log2(scmin),log2(scmax),ressc);
logscale = log2(scale(:));
detail = abs(wt.');
logpart=zeros(ressc,qres);
logFq=zeros(ressc,qres);
for nq = 1:length(q)
    for ns=1:length(s);
        DetPowQ = detail(:,ns).^q(nq) ;
        logpart(ns,nq) = log2(mean(DetPowQ));
        logFq(ns,nq) = log2((mean(DetPowQ))^(1/q(nq)));
        clear DetPowQ 
    end
end
tau=zeros(1,length(q));
for nq = 1:length(q) 
    slope = polyfit(logscale,logpart(:,nq),1) ;
    tau(nq) = slope(1);
end
Hq=(tau./q);
hh=diff(tau)./(q(2)-q(1));
Dh=((q(1:(end-1)).*hh)-tau(1:(end-1)))+1;
h=hh-0.7;