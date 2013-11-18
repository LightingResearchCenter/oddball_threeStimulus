%clear all,close all

function [h_loc,maxima_idx,h,h_row,h_pct,h_intr]=XWTMM_local(xsignal,ysignal,scmin,scmax,fig)

%load Signal7;
%signal=signal7(1:1024);
%N=length(signal);

%PARAMETERS-----------------------------

%scmin=2;
%scmax=N/5;
%scres=100;

%---------------------------------------

N=length(xsignal);
q=[0,2];
J1=round(log2(scmax/scmin)/(1/16));
[xCoefs,period,scale] = wavelet(xsignal,1,1,1/16,scmin,J1,'morlet',6);
[yCoefs] = wavelet(ysignal,1,1,1/16,scmin,J1,'morlet',6);
aCoefs=xCoefs.*conj(yCoefs);
Coefs=sqrt(abs(aCoefs));

logscale=log2(scale);
maxmap=zeros(size(Coefs));
maxmap_pks=zeros(size(Coefs));
partfunc=zeros(2,1);
M_pks=zeros(1,size(maxmap,1));

for n=1:size(maxmap,1);
    [pks,locs] = findpeaks(Coefs(n,:),'minpeakdistance',ceil(scale(n)/4));
    maxmap(n,locs)=1;
    maxmap_pks(n,locs)=log2(abs(pks));
    for nq=1:length(q),
        partfunc(nq)=sum(abs(pks).^q(nq));
    end
    M_pks(n)=sqrt(partfunc(2)/partfunc(1));
end

fit_idx=find(scale>2, 1 ):find(scale<50, 1, 'last' );
P=polyfit(logscale(fit_idx),log2(M_pks(fit_idx)),1);
h_fit=polyval(P,logscale);
h=zeros(size(maxmap));
h2=zeros(size(maxmap));
maxima_idx=find(maxmap(1,:)==1);
h_mean=P(1)*log2(N)+P(2);

for n=1:size(maxmap,1);
    map_idx=find(maxmap(n,:)==1);
    h2(n,map_idx)=(maxmap_pks(n,map_idx)-h_mean)./(logscale(n)-log2(N));
    h(n,map_idx)=(maxmap_pks(n,map_idx)-h_fit(n))./(logscale(n)-log2(N));
    clear map_idx
end

h_cone=zeros(size(maxmap,1),1);
h_cone_idx=zeros(size(maxmap,1),1);
h_loc=zeros(length(maxima_idx),1);

h2_row=h2(:);
h2_row(h2_row==0)=[];
h_row=h(:);
h_row(h_row==0)=[];
h_row=h_row-median(h_row)+median(h2_row);

for nn=1:N;%length(maxima_idx)
    for n=1:size(maxmap,1);
        ind_start=max([nn-floor(scale(n)/2) 1]);
        ind_stop=min([nn+floor(scale(n)/2) size(maxmap,2)]);
        cone_idx=ind_start:ind_stop;
        h_in_cone_idx=find(maxmap(n,ind_start:ind_stop)==1);
        if ~isempty(h_in_cone_idx)
           h_cone(n)=mean(h(n,cone_idx(h_in_cone_idx)));
           h_cone_idx(n)=n;
        end
        clear cone_idx
    end
    h_cone_idx(h_cone_idx==0)=[];
    h_loc(nn)=median(h_cone(h_cone_idx));
end
h_loc=h_loc-mean(h_loc)+median(h2_row);

h_pct=prctile(h_row,[2.5 25 50 75 97.5]);
h_intr=diff(h_pct);
%h_loc_pct=prctile(h_loc,[ 0.2 25 50 75 97.5]);
%h_loc_intr=diff(h_loc_pct);
%fig=1;
if fig==1;
    figure;plot(logscale,log2(M_pks),'bo',logscale,h_fit,'r-','LineWidth',2);
    title('Average scaling of the maxima lines'),
    ylabel('log(Amp)'),xlabel('log(scale)')
    
    figure;
    subplot(211)
    plot(1:N,xsignal,'b-',1:N,ysignal,'r-');
    title('Response series'); xlabel('Response number'); ylabel('Response time (ms)');
    subplot(212)
    plot(1:N,h_loc','b.',maxima_idx,h(1,maxima_idx),'r.');
    title('Local hölder regularity'); xlabel('Response number'); ylabel('h(t)');legend('median all scales', 'smallest scale')
    
    figure;
    viewWTLM(maxmap,scale,Coefs);
    title('Maxima lines for the h(t) computation');xlabel('Response number');ylabel('scale');
    
    [Pfreq,Bin]=hist(h_row,ceil(sqrt(length(h_row))));
    
    figure;
    bar(Bin,Pfreq,'r');hold all,
    plot(h_pct(1).*ones(2,1),[0 max(Pfreq)+1],'b--','LineWidth',2);hold all
    plot(h_pct(3).*ones(2,1),[0 max(Pfreq)+1],'b--','LineWidth',2);hold all
    title('Histogram of local hôlder exponent'); ylabel('P(h)'); xlabel('h');
end  