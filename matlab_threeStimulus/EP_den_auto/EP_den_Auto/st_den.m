function [YDEN,y,den_coeff]=st_den(x, den_coeff, sc, samples, handles)

denoised_coeff={};
for i=1:sc
    denoised_coeff{i} = downsample(den_coeff(i,:),2^i);
end
denoised_coeff{sc+1} = downsample(den_coeff(sc+1,:),2^sc);

[h,g,rh,rg] = load_filter_coefficients;
lx=length(x);
lh=length(h);
lg=length(g);
lrh=length(rh);
lrg=length(rg);
sweeps =length(x)/samples;
YDEN=[];
m=1; 
n=samples;

for sw=1:sweeps
    
    y=[];
    X=x(m:n);
    
    % convolution
    th=X;    
    for s=1:sc
        
        tga=zeros(1,length(th));tha=zeros(1,length(th));
        
        for k=1:length(th); % filtering
            for i=-lh/2:lh/2-1;
                if i>-k;
                    if i<length(th)-k;
                        tha(k) = tha(k) + th(i+k) * h(i+lh/2+1);
                    end; 
                end
            end 
            for i=-lg/2:lg/2-1;
                if i>-k;
                    if i<length(th)-k;
                        tga(k) = tga(k) + th(i+k) * g(i+lg/2+1);
                    end; 
                end
            end 
        end
        tg=tga; th=tha;
        tg=tg(1:2:length(tg)); % decimation
        th=th(1:2:length(th));
  
        % denoising details
        dc=denoised_coeff{s};
        tg(dc==0)=0;    

        %reconstruction and fill with 0s 
        for j=s:-1:1        
            
            if j==s
                tg=[tg;zeros(1,length(tg))];
                tg=tg(:)';

                tm=tg;
                tma=zeros(1,length(tm));
                for k=1:length(tm);
                    for i=-lrg/2+1:lrg/2;
                        if i>-k;
                            if i<length(tm)-k;
                                tma(k) = tma(k) + tm(i+k) * rg(i+lrg/2);
                            end; 
                        end
                    end 
                end
                tm = tma;
            
            else
                tg=[tg;zeros(1,length(tg))];
                tg=tg(:)';
                tm=[tm;zeros(1,length(tm))];
                tm=tm(:)';
                tma=zeros(1,length(tm));
                for k=1:length(tm);
                    for i=-lrh/2+1:lrh/2;
                        if i>-k;
                            if i<length(tm)-k;
                                tma(k) = tma(k) + tm(i+k) * rh(i+lrh/2);
                            end; 
                        end
                    end 
                end
                tm = tma;
            end
            
        end
        y(s,:)=tm;
        den_coeff(s,:)=tg;
    end

    %************************************************
    % denoising the last approximation
    dc=denoised_coeff{sc+1};
    th(dc==0)=0;

    %************************************************
    % reconstruction of the last approximation
    ta=th;
    for j=1:sc                                  
        th=[th;zeros(1,length(th))];
        th=th(:)';
        %ups_org_coeff(sc+1,:)=th;
        ta=[ta;zeros(1,length(ta))];
        ta=ta(:)';
        taa=zeros(1,length(ta));
        for k=1:length(ta);
            for i=-lrh/2+1:lrh/2;
                if i>-k;
                    if i<length(ta)-k;
                        taa(k) = taa(k) + ta(i+k) * rh(i+lrh/2);
                    end; 
                end
            end 
        end
        ta = taa;
    end

    den_coeff(sc+1,:)=th;
    y(sc+1,:)=ta;
    YDEN(sw,:)=sum(y);
    m=n+1;
    n=n+samples;
    
end   