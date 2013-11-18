function [s,q,Hq,h,Dh]=MFXMODWT(xsignal,ysignal,scmax,qmin,qmax,qres)

% Multifractal Maximum Overlap Discrete Wavelet Transformation for covariation between two variables(MFXMODWT)
%
% [s,q,Hq,h,Dh]=MFXMODWT(xsignal,ysignal,scmax,qmin,qmax,qres)
%
% INPUT PARAMETERS---------------------------------------------
% 
% xsignal:      input signal 1
% ysignal:      input signal 2
% scmax:        upper bound of the window size s
% qmin          lower bound of q
% qmax          upper bound of q
% qres          number of elements in q
%
% OUTPUT VARIABLES--------------------------------------------
% 
% s:            scale
% q:            q-order 
% Hq:           q-generalized Hurst exponent
% h:            Hölder exponent 
% Dh:           Multifractal spectrum 
%
% EXAMPLE-----------------------------------------------------
%
% T=40.96 ; rmin=0.02 ; Dt=rmin/2 ; Law=0 ; ParamLaw1=0.2; ParamLaw2=0; C=1 ; q1=(-10:10)/2;
% [Qr1,time,TtheoQ,q] = IDCnoise_epl(T,rmin,Dt,Law,ParamLaw1,C,q1); 
% H=0.75 ; oversamprate = 16 ; printout = 0;
% [Ar,VH1,BH] = synthArVH(H,Qr1,Dt,oversamprate,printout) ;
% [Ar,VH2,BH] = synthArVH(H,Qr1,Dt,oversamprate,printout) ;
% [s,q,Hq,h,Dh]=MFXMODWT(diff(VH1),diff(VH2),32,0.1,3,30);
% figure;
% subplot(211)
% plot(1:4096,diff(VH1),'b-',1:4096,diff(VH2),'r-');xlabel('time');
% ylabel('amplitude'); title('multiplicative cascading noise \DeltaB_H(A(t))and fractional Gaussian noise \DeltaB_H(t)');
% legend('\DeltaB_H(A(t))','\DeltaB_H(t)');
% subplot(223)
% plot(q,Hq1,'bo-',q,Hq2,'ro-');xlabel('q'); ylabel('H(q)'); title('q-generalized Hurst exponent');
% legend('\DeltaB_H(A(t))','\DeltaB_H(t)');
% subplot(224)
% plot(h1, Dh1,'bo-',h2, Dh2,'ro-');xlabel('h(q)');ylabel('D(h)');title('multifractal half-spectrum');
% legend('\DeltaB_H(A(t))','\DeltaB_H(t)');
%
% ---------------------------------------------------------------
% Written by Espen A. F. Ihlen (espen.ihlen@ntnu.no), 2010

if length(xsignal)~=length(ysignal),
    disp('---ERROR: xsignal and ysignal MUST be the same length--');
end

log_scmin=1;
log_scmax=floor(log2(scmax));
q=linspace(qmin,qmax,qres);

s=2.^(log_scmin:log_scmax);
loge = log2(exp(1)) ;
q = sort(q);
wtf='la20';
xsignal=xsignal-mean(xsignal);
ysignal=ysignal-mean(ysignal);

[xWJt, xVJt, xatt]=modwt(xsignal,wtf,log_scmax,'reflection',[]);
[xDJt,xSJt] = imodwt_mra(xWJt,xVJt,xatt);

[yWJt, yVJt, yatt]=modwt(ysignal,wtf,log_scmax,'reflection',[]);
[yDJt,ySJt] = imodwt_mra(yWJt,yVJt,yatt);

wave_coeff=zeros(length(xsignal),length(s));
for ns=1:length(s),
    X=hilbert(xDJt(1:length(xsignal),ns));
    Y=hilbert(yDJt(1:length(ysignal),ns));
    wave_coeff(1:length(xsignal),ns)=abs(X.*conj(Y));
    clear X Y
end
    
for j=1:log_scmax
    nj(j) = length(wave_coeff);
    for  k=1:length(q)    
        q_wave_coeff = wave_coeff(:,j).^(q(k)/2) ;            % Eq. (A13)
        logpart = log2(mean(q_wave_coeff));                   % Eq. (A13)   
        V_q = std(q_wave_coeff)^2/mean(q_wave_coeff)^2 ;      % Eq. (A17)      
        Var1(k,j) =  loge^2 * V_q /nj(j) ;                    % X_q in Eq. (A16)  
        gj = -loge/2 * V_q /nj(j) ;                      
        Var2(k,j) = logpart - gj ;   
        Var2(k,j) = Var2(k,j) - j*q(k)/2;                     % Y_q in Eq. (A16)
    end
end
for k=1:length(q);
    log_scale = log_scmin:log_scmax ;
    var_x(log_scale) = Var1(k,log_scale)';
    X_q = var_x; 
    njj  = nj(log_scale);
    var_y(log_scale)  = Var2(k,log_scale)';
    Y_q=var_y;
    
    W1 = sum(1./X_q) ;                 %weight 1 in Eq. (A15)
    W2 = sum(log_scale./X_q) ;         %weight 2 in Eq. (A15)
    W3 = sum(log_scale.^2./X_q) ;      %weight 3 in Eq. (A15)
    slope(k)  = sum(((W1 * log_scale - W2) ./ X_q / (W1*W3-W2*W2)) .* Y_q); %Eq. (A14)                   
end

tau=slope;
Hq=tau./q+1.5;
hh=diff(tau)./(q(2)-q(1));
Dh=((q(1:(end-1)).*hh)-tau(1:(end-1)))+1;
h=hh+0.5;
