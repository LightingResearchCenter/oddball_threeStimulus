function [s,q,Hq,h,Dh,logFq]=MFDXA_modwt(xsignal,ysignal,qmin,qmax,qres)
%
% Multifractal detrended crosscorrelation analysis (MFDFA) with MODWT based detrending
%
% [s,q,Hq,h,Dh,logFq]=MFDXA_modwt(xsignal,ysignal,qmin,qmax,qres);
%
% INPUT PARAMETERS-----------------------------------------------------
%
% xsignal:      input signal 1
% ysignal:      input signal 2
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
% ---------------------------------------------------------------
% Written by Espen A. F. Ihlen (espen.ihlen@ntnu.no), 2010

if length(xsignal)~=length(ysignal),
    disp('---ERROR: xsignal and ysignal is not the same length--');
end
warning off
xsignal=cumsum(xsignal-mean(xsignal));
ysignal=cumsum(ysignal-mean(ysignal));
N=length(xsignal);

[xDJt,xSJt]=modwt_mra(xsignal,'la20', floor(log2(length(xsignal))),'reflection');
[yDJt,ySJt]=modwt_mra(ysignal,'la20', floor(log2(length(ysignal))),'reflection');
xFluct=zeros(size(xDJt,2),length(xsignal));
yFluct=zeros(size(yDJt,2),length(ysignal));

for n=1:size(xDJt,2)
    xFluct(n,:)=xsignal'-(sum(xDJt(1:N,n+1:end),2)+xSJt(1:N));
    yFluct(n,:)=ysignal'-(sum(yDJt(1:N,n+1:end),2)+ySJt(1:N));
end
xFluctRev=fliplr(xFluct);
yFluctRev=fliplr(yFluct);

xan=zeros(size(xDJt,2),length(xFluct));
xphi=zeros(size(xDJt,2),length(xFluct));
xinstfreq=zeros(size(xDJt,2),length(xFluct));

yan=zeros(size(yDJt,2),length(yFluct));
yphi=zeros(size(yDJt,2),length(yFluct));
yinstfreq=zeros(size(yDJt,2),length(yFluct));

deltat=1;

for n=1:size(xDJt,2);
    xan(n,:)=hilbert(xDJt(1:N,n)')';
    xphi(n,:)=angle(xan(n,:));
    xphi(n,:)=unwrap(xphi(n,:));
    xinstfreq(n,:)=phasediff(xphi(n,:),deltat)./(2*pi);
end
for n=1:size(yDJt,2);
    yan(n,:)=hilbert(yDJt(1:N,n)')';
    yphi(n,:)=angle(yan(n,:));
    yphi(n,:)=unwrap(yphi(n,:));
    yinstfreq(n,:)=phasediff(yphi(n,:),deltat)./(2*pi);
end

s=round(1./abs(median((yinstfreq+xinstfreq)./2,2)));
s(find(s>length(xsignal)))=[];

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
        xSeg=xFluct(ns,SegNumb);                                             
        xSegRev=xFluctRev(ns,SegNumb);  
        ySeg=yFluct(ns,SegNumb);                                             
        ySegRev=yFluctRev(ns,SegNumb);
        for nq=1:length(q);
            Var(v,nq)=sum(((xSeg.^2).*(ySeg.^2))/s(ns))^(q(nq)/4);                       
            Varr(v,nq)=sum(((xSegRev.^2).*(ySegRev.^2))/s(ns))^(q(nq)/4); 
        end;
        clear SegNumb xSeg xSegRev ySeg ySegRev 
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
           