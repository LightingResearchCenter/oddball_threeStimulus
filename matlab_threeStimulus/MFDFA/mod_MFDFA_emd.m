function [s,q,Hq,h,Dh,logFq,Sdt,C]=mod_MFDFA_emd(signal,qmin,qmax,qres,qbase)
%
% Non-scale invariant version of Multifractal detrended fluctuation analysis (MFDFA) with EMD based
% detrending
%
% [s,q,Hq,h,Dh,logFq,Sdt,C]=mod_MFDFA_emd(signal,qmin,qmax,qres,qbase);
%
% INPUT PARAMETERS-----------------------------------------------------
%
% signal:       input signal
% qmin          lower bound of q
% qmax          upper bound of q
% qres          number of elements in q
% qbase         baseline q use to compute Sdt and C (q = 0 recommended)
%
% OUTPUT VARIABLES-----------------------------------------------------
%
% s:            scale
% q:            q-order 
% Hq:           q-generalized Hurst exponent
% h:            Hölder exponent 
% Dh:           Multifractal spectrum 
% logFq:        q-generalized scaling function F(s,q) in log coordinates 
% Sdt           q-average scale dependency for logFq
% C             the q dependent intercepts
%
% ---------------------------------------------------------------
% Written by Espen A. F. Ihlen (espen.ihlen@ntnu.no), 2010

warning('off');
signal=cumsum(signal-mean(signal));
N=length(signal);
imf=emd(signal);
Fluct=zeros(size(imf,1)-1,length(signal));
for n=1:size(imf,1)-1
    Fluct(n,:)=signal-sum(imf(n+1:end,:),1);
end
FluctRev=fliplr(Fluct);

an=zeros(size(imf,1)-1,length(Fluct));
phi=zeros(size(imf,1)-1,length(Fluct));
instfreq=zeros(size(imf,1)-1,length(Fluct));

deltat=1;

for n=1:size(imf,1)-1;
    an(n,:)=hilbert(imf(n,:)')';
    phi(n,:)=angle(an(n,:));
    phi(n,:)=unwrap(phi(n,:));
    instfreq(n,:)=phasediff(phi(n,:),deltat)./(2*pi);
end

s=round(1./abs(mean(instfreq,2)));
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
qnumb=find(q==qbase);
C=zeros(1,length(q));
%Arg=zeros(ressc,qres);
Fq_mod=zeros(2,qres);

for nq=1:length(q);
    Fq_mod(:,nq)=polyfit(log2(Fq(:,qnumb)),log2(Fq(:,nq)),1);
    %Hq_mod(nq)=Fq_mod(1,nq);
    P=polyfit(log2(s),log2(Fq(:,nq)),1);
    %Hq(nq)=P(1);
    if q(nq)==qbase;
       Kp=P(2);
    end
end;

for nq=1:length(q);
    Hq(nq)=Fq_mod(1,nq);
    C(nq)=mean(log2(Fq(:,nq))-(Fq_mod(1,nq).*log2(Fq(:,qnumb))));
    for ns=1:length(s);
        Arg(ns,nq)=(1/Fq_mod(1,nq))*(log2(Fq(ns,nq))-C(nq));
    end;
end;
Sdt(1:length(s))=((1/Hq(qnumb)).*mean(Arg,2))+Kp;

tau=(q.*Hq)-1;

h=diff(tau)./(q(2)-q(1));
Dh=(q(1:(end-1)).*h)-tau(1:(end-1));

           
           