function [s,q,Hq,h,Dh,logFq,Sdt,C]=mod_MFDFA(signal,m,scmin,scmax,ressc,qmin,qmax,qres,qbase)

% non-scale inveriant Multifractal detrended fluctuation analysis (MFDFA)
%
% [s,q,Hq_mod,Fq,alpha_mod,f_alpha_mod,Scaling,Coeff]=mod_MFDFA(signal,m,scmin,scmax,ressc,qmin,qmax,qres,qbase);
%
% INPUT PARAMETERS---------------------------------------------
%
% signal:       input signal
% m:            polynomial order for the detrending
% scmin:        lower bound of the window size s
% scmax:        upper bound of the window size s
% ressc:        number of elements in s 
% qmin          lower bound of q
% qmax          upper bound of q
% qres          number of elements in q
% qbase         referrance q (as small as possible)
%
% OUTPUT PARAMETERS--------------------------------------------
%
% s:            scale
% q:            q-order 
% Hq_mod:       q-generalized Hurst exponent (non-scale inveriant)
% h_mod:        hölder exponent (non-scale invariant)
% Dh_mod:       Multifractal spectrum (non_scale invariant)
% logFq:        q-generalized scaling function Fq(s,q) in log2 coordinates
% Scaling:      The non-invariant scaling (q-average)
% Coeff:        q-dependent intersections
%
% --------------------------------------------------------------
% 
% Written by Espen A. F. Ihlen (espenale@svt.ntnu.no), 2009

Fluct=cumsum(signal-mean(signal)./std(signal));
FluctRev=fliplr(Fluct);
N=length(Fluct);

ScaleNumb=linspace(log2(scmin),log2(scmax),ressc);
s=round(2.^ScaleNumb);
q=linspace(qmin,qmax,qres);
znumb=find(q==0);
Fq=zeros(length(s),length(q));
for ns=1:length(s),                                                      
    Ns=floor(s(ns)\N);
    Var=zeros(Ns,length(q));
    Varr=zeros(Ns,length(q));
    for v=1:Ns,                                                         
        SegNumb=((((v-1)*s(ns))+1):(v*s(ns)))';
        Seg=Fluct(SegNumb);                                             
        SegRev=FluctRev(SegNumb);                                      
        poly=polyfit(SegNumb,Seg',m);                                   
        polyr=polyfit(SegNumb,SegRev',m);                               
        fit=polyval(poly,SegNumb);                                     
        fitr=polyval(polyr,SegNumb);        
        for nq=1:length(q);
            Var(v,nq)=((sum((Seg'-fit).^2))/s(ns))^(q(nq)/2);                       
            Varr(v,nq)=((sum((SegRev'-fitr).^2))/s(ns))^(q(nq)/2); 
        end;
        clear SegNumb Seg SedRev poly polyr fit fitr
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
Arg=zeros(ressc,qres);
Fq_mod=zeros(2,qres);

for nq=1:length(q);
    Fq_mod(:,nq)=polyfit(log2(Fq(:,qnumb)),log2(Fq(:,nq)),1);
    %Hq_mod(nq)=Fq_mod(1,nq);
    P=polyfit(log2(s'),log2(Fq(:,nq)),1);
    Hq(nq)=P(1);
    if q(nq)==qbase;
       Kp=P(2);
    end
end;

%Hq_mod=Hq_mod+Hq(qnumb)-1;

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

                      