function [s,q,Hq,h,Dh,logFq]=MFDFA(signal,m,scmin,scmax,ressc,qmin,qmax,qres)

%{
    debugMatFileName = 'MFDFA.mat';
    if nargin == 0
        load(debugMatFileName)
    else
        save(debugMatFileName)        
    end  
%}


%
% Multifractal detrended fluctuation analysis (MFDFA)
%
% [s,q,Hq,h,Dh,logFq]=MFDFA(signal,m,scmin,scmax,ressc,qmin,qmax,qres);
%
% INPUT PARAMETERS-----------------------------------------------------
%
% signal:       input signal
% m:            polynomial order for the detrending
% scmin:        lower bound of the window size s
% scmax:        upper bound of the window size s
% ressc:        number of elements in s 
% qmin          lower bound of q
% qmax          upper bound of q
% qres          number of elements in q
%
% OUTPUT VARIABLES-----------------------------------------------------
%
% s:            scale
% q:            q-order 
% Hq:           q-generalized Hurst exponent
% logFq:        q-generalized scaling function F(s,q) in log coordinates 
% lin_fit       linear least square fit to logFq
% h:            H�lder exponent 
% Dh:           Multifractal spectrum 
%
% EXAMPLE-----------------------------------------------------
%
% T=40.96 ; rmin=0.02 ; Dt=rmin/2 ; Law=0 ; ParamLaw1=0.2; ParamLaw2=0; C=1 ; q1=(-10:10)/2;
% [Qr1,time,TtheoQ,q] = IDCnoise_epl(T,rmin,Dt,Law,ParamLaw1,C,q1); 
% [Qr2,time,TtheoQ,q] = IDCnoise_epl(T,rmin,Dt,Law,ParamLaw2,C,q1); 
% H=0.75 ; oversamprate = 16 ; printout = 0;
% [Ar,VH1,BH] = synthArVH(H,Qr1,Dt,oversamprate,printout) ;
% [Ar,VH2,BH] = synthArVH(H,Qr2,Dt,oversamprate,printout) ;
% [s,q,Hq1,h1,Dh1,logFq1]=MFDFA(diff(VH1),1,10,200,40,-3,3,61);
% [s,q,Hq2,h2,Dh2,logFq2]=MFDFA(diff(VH2),1,10,200,40,-3,3,61);
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

debug_MFDFA = 0;

% debug param
if debug_MFDFA == 1
    rows = 3;
    cols = 5;
    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    fig = figure('Color', 'w');
        set(fig, 'Position', [0.05*scrsz(3) 0.20*scrsz(4) 0.92*scrsz(3) 0.62*scrsz(4)])
end

%warning off;
Fluct=cumsum(signal-mean(signal)./std(signal));
FluctRev=fliplr(Fluct);
N=length(Fluct);

ScaleNumb=linspace(log2(scmin),log2(scmax),ressc);
s=round(2.^ScaleNumb);
q=linspace(qmin,qmax,qres);
znumb=find(q==0);
Fq=zeros(length(s),length(q));

if debug_MFDFA == 1
    subplot(rows,cols,1); plot(Fluct); title('Fluct'); % debug plot
    subplot(rows,cols,2); plot(s); title('s'); % debug plot
    drawnow
end

for ns=1:length(s),%disp(strcat('computing scale number_',num2str(ns)));                                                      
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
        
        if debug_MFDFA == 1
            subplot(rows,cols,3); plot(Var); title('Var'); % debug plot 
            subplot(rows,cols,4); plot(Varr); title('Varr'); % debug plot
            drawnow        
        end
        clear SegNumb Seg SedRev poly polyr fit fitr
        
    end
    
    for nq=1:length(q),
        Fq(ns,nq)=((sum(Var(:,nq))+sum(Varr(:,nq)))/(2*Ns))^(1/q(nq));
    end
    
    if debug_MFDFA == 1
        subplot(rows,cols,5); plot(Fq); title('Fq 1'); % debug plot
    end
    Fq(ns,znumb)=(Fq(ns,znumb-1)+Fq(ns,znumb+1))./2;

    if debug_MFDFA == 1
        subplot(rows,cols,6); plot(Fq); title('Fq znumb'); % debug plot
        drawnow    
    end
    clear Var Varr
end

logFq=log2(Fq);
Hq=zeros(1,length(q));
lin_fit=zeros(length(s),length(q));

if debug_MFDFA == 1
    subplot(rows,cols,7); plot(logFq); title('logFq'); % debug plot
end

for nq=1:length(q);
    
    % remove the NaN values if any out there
    s = s(~isnan(logFq(:,nq)));
    logFq2(:,nq) = logFq((~isnan(logFq(:,nq))),nq);
   
    if debug_MFDFA == 1
        subplot(rows,cols,8); plot(log2(s'), logFq2(:,nq))
        drawnow
    end
    P = polyfit(log2(s'),logFq2(:,nq),1);    
    lin_fit(1:length(s),nq)=polyval(P,log2(s'));
    Hq(nq)=P(1);    
end

if debug_MFDFA == 1
    subplot(rows,cols,9); plot(Hq(nq)); title(['Hq, ', 'noNaNs = ', num2str(sum(isnan(Hq))), ', noINFs = ', num2str(sum(isinf(Hq)))]); % debug plot
    subplot(rows,cols,10); plot(lin_fit); title('lin_fit'); % debug plot
end

tau=(q.*Hq)-1;

hh=diff(tau)./(q(2)-q(1));
Dh=(q(1:(end-1)).*hh)-tau(1:(end-1));
h=hh-1;
   
if debug_MFDFA == 1

    subplot(rows,cols,11); plot(hh); title('hh'); % debug plot
    subplot(rows,cols,12); plot(Dh); title('Dh'); % debug plot
    subplot(rows,cols,13); plot(h); title('h'); % debug plot

    fileNameOut = sprintf('%s%s', 'debug_MFDFA_steps2', '.png');
    export_fig(fileNameOut, '-r400', '-a1', fig)
end