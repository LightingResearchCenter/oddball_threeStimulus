function [coeff,denav,den_coeff,y,yo]=Run_Neigh(av,handles);
stim=handles.par.stim;
sc=handles.par.scales;
[h,g,rh,rg] = load_filter_coefficients;
lav=length(av);lh=length(h);lg=length(g);lrh=length(rh);lrg=length(rg);
y=zeros(sc+1,lav);
yo=zeros(sc+1,lav);
coeff=zeros(sc+1,lav);
den_coeff=zeros(sc+1,lav);
%% CONVOLUTION
th=av;
for s=1:sc
    tga=zeros(1,length(th));tha=zeros(1,length(th));
    for k=1:length(th);                    %filtering
        for i=-lh/2:lh/2-1;
            if i>-k;if i<length(th)-k;
                    tha(k) = tha(k) + th(i+k) * h(i+lh/2+1);
            end; end
        end 
        for i=-lg/2:lg/2-1;
            if i>-k;if i<length(th)-k;
                    tga(k) = tga(k) + th(i+k) * g(i+lg/2+1);
            end; end
        end 
    end
    tg=tga; th=tha;
    tg=tg(1:2:length(tg));                 %decimation
    th=th(1:2:length(th));
    TG=tg;
    TH=th;
%% denoising details
    L_tg=length(tg);
%     tgm=tg(1:L_tg/2)-mean(tg(1:L_tg/2));  
    Pow=nextpow2(stim-1);
    L_baseline=2^Pow;
    tgm=tg(1:L_baseline/(2^s))-mean(tg(1:L_baseline/(2^s)));
    bl_var=median(abs(tgm)/0.6745);   % estimate basline noise
    LS = [tg(2:length(tg)) tg(1)];
	RS= [tg(length(tg)) tg(1: length(tg)-1)];
    e_tg=RS.^2+tg.^2+LS.^2;
    thresh(s) =(sqrt(2*(bl_var^2)* log(L_tg)))^2; 
    e_tg( find( abs(e_tg)<=abs(thresh(s)) ) ) = 0;
    mask=e_tg & ones(1,L_tg);
    tg=mask.*tg;
%% reconstruction and fill with 0s 
    for j=s:-1:1                          
        if j==s
            TG=[TG;zeros(1,length(TG))];
            TG=TG(:)';
            tg=[tg;zeros(1,length(tg))];
            tg=tg(:)';
            
            tm=tg;                                       % reconstruct the denoised signal
            tma=zeros(1,length(tm));  
            for k=1:length(tm);
                for i=-lrg/2+1:lrg/2;
                    if i>-k;if i<length(tm)-k;
                            tma(k) = tma(k) + tm(i+k) * rg(i+lrg/2);
                    end; end
                end
            end
            tm = tma;
            tmo=TG;                                      % reconstruct the bands of the original signal
            tma=zeros(1,length(tmo));
            for k=1:length(tmo);
                for i=-lrg/2+1:lrg/2;
                    if i>-k;if i<length(tmo)-k;
                            tma(k) = tma(k) + tmo(i+k) * rg(i+lrg/2);
                    end; end
                end
            end
            tmo=tma;
        else
            TG=[TG;zeros(1,length(TG))];
            TG=TG(:)';
            tg=[tg;zeros(1,length(tg))];
            tg=tg(:)';
            
            tm=[tm;zeros(1,length(tm))];
            tm=tm(:)';
            tma=zeros(1,length(tm));
            for k=1:length(tm);
                for i=-lrh/2+1:lrh/2;
                    if i>-k;if i<length(tm)-k;
                            tma(k) = tma(k) + tm(i+k) * rh(i+lrh/2);
                    end; end
                end 
            end
            tm = tma;
            tmo=[tmo;zeros(1,length(tmo))];
            tmo=tmo(:)';
            tma=zeros(1,length(tmo));
            for k=1:length(tmo);
                for i=-lrh/2+1:lrh/2;
                    if i>-k;if i<length(tmo)-k;
                            tma(k) = tma(k) + tmo(i+k) * rh(i+lrh/2);
                    end; end
                end 
            end
            tmo=tma;
            end
    end
   y(s,:)=tm;
   yo(s,:)=tmo;
   coeff(s,:)=TG;
   den_coeff(s,:)=tg;
end

%% denoising the last approximation
    L_th=length(th);
%     thm=th(1:L_th/2)-mean(th(1:L_th)/2);
    thm=th(1:L_baseline/(2^s))-mean(th(1:L_baseline/(2^s)));
    bl_var=median(abs(thm)/0.6745);
    LS = [th(2:length(th)) th(1)];
	RS= [th(length(th)) th(1: length(th)-1)];
    e_th=RS.^2+th.^2+LS.^2;     
    thresh(sc+1) =(sqrt(2*(bl_var^2)* log(L_th)))^2;   
    e_th( find( abs(e_th)<=abs(thresh(sc+1)) ) ) = 0;
    mask=e_th & ones(1,L_th);
    th=mask.*th;
%% reconstruction of the last approximation
    ta=th;
    tao=TH;
    for j=1:sc 
        TH=[TH;zeros(1,length(TH))];
        TH=TH(:)';
        th=[th;zeros(1,length(th))];
        th=th(:)';
        ta=[ta;zeros(1,length(ta))];
        ta=ta(:)';
        taa=zeros(1,length(ta));
        for k=1:length(ta);
            for i=-lrh/2+1:lrh/2;
                   if i>-k;if i<length(ta)-k;
                             taa(k) = taa(k) + ta(i+k) * rh(i+lrh/2);
                   end; end
            end
        end
        ta = taa;
        tao=[tao;zeros(1,length(tao))];
        tao=tao(:)';
        taa=zeros(1,length(tao));
        for k=1:length(tao);
            for i=-lrh/2+1:lrh/2;
                   if i>-k;if i<length(tao)-k;
                             taa(k) = taa(k) + tao(i+k) * rh(i+lrh/2);
                   end; end
            end
        end
    tao=taa;    
    end
y(sc+1,:)=ta;
yo(sc+1,:)=tao;
coeff(sc+1,:)=TH;
den_coeff(sc+1,:)=th;
denav=sum(y);
