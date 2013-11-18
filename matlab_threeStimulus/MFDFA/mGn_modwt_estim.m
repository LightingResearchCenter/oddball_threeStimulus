function [s,time,Ht,Betat,logFt]=mGn_modwt_estim(signal,win,scmax,wtf)

% MODWT based estimation of H(t) 
%
% [s,time,Ht,Betat,logFt]=mGn_modwt_estim(signal,win,scmax,wtf)
%
% INPUT ARGUMENTS---------------------------------------------
%
% signal:       input signal
% win:          sample size of the moving time window where H(t)is computed
% scmax:        upper bound of the window size s
% wtf:          waveform (e.g., 'la4','la8','la16');
%
% OUTPUT VARIABLES--------------------------------------------
%
% s:            scale
% time:         time indices
% Ht:           local Hurst exponent H(t)
% Betat:        local Beta exponent
% logFt:        local scaling function F(s,t)
%
% EXAMPLE-----------------------------------------------------
% 
% N=1024;
% t=1:N;
% Ht=0.5+0.5.*(sin(0.0025.*pi.*t));
% [mBm,mGn]=mBm_mGn(N,Ht);
% [s,time,Ht_est,Betat,logFt]=mGn_modwt_estim(mGn,150,32,'la8');
% figure;
% plot(time,Ht_est,'b-',time,Ht(time),'r-')
% xlabel('time');ylabel('H'),legend('Estimated H(t)','H(t)')
%
%--------------------------------------------------------------
% Written by Espen A. F. Ihlen (espenale@svt.ntnu.no), 2009

J0=floor(log2(scmax));
N=length(signal);

[WJt, VJt, att]=modwt(signal,wtf,J0,'reflection',[]);
TWJt=modwt_cir_shift(WJt,[],wtf);


if win>1
    [rwvar, CI_rwvar, time] = modwt_running_wvar(TWJt(:,1:J0),[ceil(win/2)+1:N-ceil(win/2)],1,win,'chi2eta3','biased',wtf);
else
    rwvar=TWJt(1:N,:).^2;
end    
    
rwvar1=zeros(length(rwvar),J0);
CI_rwvar1=zeros(J0,2,length(rwvar));
Betat=zeros(length(rwvar),1);
Ht=zeros(length(rwvar),1);
logFt=zeros(J0,length(rwvar));

if win>1
    for n_var=1:length(rwvar),
        CI_rwvar1(1:J0,1:2,n_var)=CI_rwvar(n_var,1:J0,1:2);
        rwvar1(1:J0,n_var)=rwvar(n_var,1:J0);
        [CJr,fr,CI_CJr,f_bandr]=modwt_wvar_psd(rwvar1(1:J0,n_var),1,CI_rwvar1(1:J0,1:2,n_var));
        P=polyfit(-log2(fr),log2(CJr),1);
        Betat(n_var)=P(1);
        Ht(n_var)=(P(1)+1)/2;
        logFt(1:length(CJr),n_var)=log2(CJr);
        s=1./fr;
        clear CJr
    end
else
    for n_var=1:length(rwvar),
        rwvar1 = rwvar(n_var,:);
        rwvar1 = rwvar1(:);
        j_range = [1:length(rwvar1)]';
        f_band1 = 1 ./ (2.^(j_range+1));
        f_band2 = 1 ./ (2.^(j_range));
        fr = sqrt(f_band1 .* f_band2);
        CJr = rwvar1.* 2.^j_range;
        P=polyfit(-log2(fr),log2(CJr),1);
        Betat(n_var)=P(1);
        Ht(n_var)=(P(1)+1)/2;
        logFt(1:length(CJr),n_var)=log2(CJr);
        clear CJr
    end
    s=1./fr;
end
    
    
