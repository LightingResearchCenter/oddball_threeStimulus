function [s,time,Ht,logFt]=mGn_dfa_estim(signal,win,scmin,scmax,ressc,m)

% Estimation of H(t) by windowed detrended fluctuation analysis 
%
% [s,time,Ht,logFt]=mGn_dfa_estim(signal,win,scmin,scmax,ressc,m)
%
% INPUT ARGUMENTS---------------------------------------------
%
% signal:       input signal
% win:          sample size of the moving time window where H(t)is computed
% scmin:        lower bound of the window size s
% scmax:        upper bound of the window size s
% ressc:        scale resolution
% m:            polynomial order for the detrending
%
% OUTPUT VARIABLES--------------------------------------------
%
% s:            scale
% time:         time indices
% Ht:           local Hurst exponent H(t)
% logFt:        local scaling function F(s,t)
%
% EXAMPLE-----------------------------------------------------
%
% N=1024;
% t=1:N;
% Ht=0.5+0.5.*(sin(0.0025.*pi.*t));
% [mBm,mGn]=mBm_mGn(N,Ht);
% [s,time,Ht_est,logFt]=mGn_dfa_estim(mGn,200,10,32,12,1);
% figure;
% plot(time,Ht_est(time),'b-',time,Ht(time),'r-')
% xlabel('time');ylabel('H'),legend('Estimated H(t)','H(t)')
%
% ------------------------------------------------------------
% NB: This script has a long computation time (several minutes)
%
% Written by Espen A. F. Ihlen (espenale@svt.ntnu.no), 2009

signal=signal';
s=linspace(scmin,scmax,ressc);
time=(win/2)+1:length(signal)-(win/2);
Ft=zeros(length(s),length(time));
Ht=zeros(length(time),1);

for t=(win./2)+1:length(signal)-(win./2);
    Ind=t-(win/2):t+(win/2);
    Fluct=cumsum((signal(Ind)-mean(signal(Ind)))./std(signal(Ind)));
    FluctRev=fliplr(Fluct);
    N=length(Fluct);
    for ns=1:length(s),                                                      
        Ns=floor(s(ns)\N);
        Var=zeros(Ns,1);
        Varr=zeros(Ns,1);
        for v=1:Ns,                                                         
            SegNumb=((((v-1)*s(ns))+1):(v*s(ns)))';
            Seg=Fluct(SegNumb);                                             
            SegRev=FluctRev(SegNumb);                                      
            poly=polyfit(SegNumb,Seg',m);                                   
            polyr=polyfit(SegNumb,SegRev',m);                               
            fit=polyval(poly,SegNumb);                                     
            fitr=polyval(polyr,SegNumb);        
            Var(v)=(sum((Seg'-fit).^2))/s(ns);                       
            Varr(v)=(sum((SegRev'-fitr).^2))/s(ns); 

            clear SegNumb Seg SedRev poly polyr fit fitr
        end
        Ft(ns,t)=((sum(Var)+sum(Varr))/(2*Ns))^(1/2);

        clear Var Varr
    end

    P=polyfit(log2(s'),log2(Ft(:,t)),1);
    Ht(t)=P(1);
    clear Ind Fluct FluctRev N P
end    

logFt(1:length(s),1:length(time))=log2(Ft(:,time));