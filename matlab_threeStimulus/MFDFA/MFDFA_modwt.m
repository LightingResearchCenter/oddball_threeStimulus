function [s,q,Hq,h,Dh,logFq]=MFDFA_modwt(signal,qmin,qmax,qres)
%
% Multifractal detrended fluctuation analysis (MFDFA) with MODWT based
% detrending
%
% [s,q,Hq,h,Dh,logFq]=MFDFA_modwt(signal,qmin,qmax,qres);
%
% INPUT PARAMETERS-----------------------------------------------------
%
% signal:       input signal
% qmin          lower bound of q
% qmax          upper bound of q
% qres          number of elements in q
%
% OUTPUT VARIABLES-----------------------------------------------------
%
% s:            scale
% q:            q-order 
% Hq:           q-generalized Hurst exponent
% h:            Hölder exponent 
% Dh:           Multifractal spectrum 
% logFq:        q-generalized scaling function F(s,q) in log coordinates
%
% EXAMPLE-----------------------------------------------------
%
% T=40.96 ; rmin=0.02 ; Dt=rmin/2 ; Law=0 ; ParamLaw1=0.2; ParamLaw2=0; C=1 ; q1=(-10:10)/2;
% [Qr1,time,TtheoQ,q] = IDCnoise_epl(T,rmin,Dt,Law,ParamLaw1,C,q1); 
% [Qr2,time,TtheoQ,q] = IDCnoise_epl(T,rmin,Dt,Law,ParamLaw2,C,q1); 
% H=0.75 ; oversamprate = 16 ; printout = 0;
% [Ar,VH1,BH] = synthArVH(H,Qr1,Dt,oversamprate,printout) ;
% [Ar,VH2,BH] = synthArVH(H,Qr2,Dt,oversamprate,printout) ;
% [s,q,Hq1,h1,Dh1,logFq1]=MFDFA_modwt(diff(VH1),-3,3,31);
% [s,q,Hq2,h2,Dh2,logFq2]=MFDFA_modwt(diff(VH2),-3,3,31);
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
% Written by Espen A. F. Ihlen (espen.ihlen@ntnu.no), 2010

warning('off');
signal=cumsum(signal-mean(signal));
N=length(signal);

[DJt,SJt]=modwt_mra(signal,'la20', floor(log2(length(signal))),'reflection');
Fluct=zeros(size(DJt,2),length(signal));
for n=1:size(DJt,2)
    Fluct(n,:)=signal'-(sum(DJt(1:N,n+1:end),2)+SJt(1:N));
end
FluctRev=fliplr(Fluct);

an=zeros(size(DJt,2),length(Fluct));
phi=zeros(size(DJt,2),length(Fluct));
instfreq=zeros(size(DJt,2),length(Fluct));

deltat=1;

for n=1:size(DJt,2);
    an(n,:)=hilbert(DJt(1:N,n)')';
    phi(n,:)=angle(an(n,:));
    phi(n,:)=unwrap(phi(n,:));
    instfreq(n,:)=phasediff(phi(n,:),deltat)./(2*pi);
end

s=round(1./abs(median(instfreq,2)));
s(find(s>length(signal)))=[];

q=linspace(qmin,qmax,qres);
znumb=find(q==0);
Fq=zeros(length(s),length(q));
for ns=1:length(s),%disp(strcat('computing scale number_',num2str(ns)));                                                      
    if s(ns)>N/2;
        Ns=1;
    else
        Ns=floor(s(ns)\N);
    end
    Var=zeros(Ns,length(q));
    Varr=zeros(Ns,length(q));
    for v=1:Ns,                                                         
        SegNumb=((((v-1)*s(ns))+1):(v*s(ns)))';
        Seg=Fluct(ns,SegNumb);                                             
        SegRev=FluctRev(ns,SegNumb);                                                                         
        for nq=1:length(q);
            Var(v,nq)=((sum(Seg.^2))/s(ns))^(q(nq)/2);                       
            Varr(v,nq)=((sum(SegRev.^2))/s(ns))^(q(nq)/2); 
        end;
        clear SegNumb Seg SegRev
    end
    for nq=1:length(q),
        Fq(ns,nq)=((sum(Var(:,nq))+sum(Varr(:,nq)))/(2*Ns))^(1/q(nq));
    end
    Fq(ns,znumb)=(Fq(ns,znumb-1)+Fq(ns,znumb+1))./2;
    clear Var Varr
end
logFq=log2(Fq);
Hq=zeros(1,length(q));
lin_fit=zeros(length(s),length(q));
for nq=1:length(q);
    P=polyfit(log2(s),logFq(:,nq),1);
    lin_fit(1:length(s),nq)=polyval(P,log2(s'));
    Hq(nq)=P(1);
end;

tau=(q.*Hq)-1;

hh=diff(tau)./(q(2)-q(1));
Dh=(q(1:(end-1)).*hh)-tau(1:(end-1));
h=hh-1;
           