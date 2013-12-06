
function [EMDdenoised, IMF] = EMDdenoise_modPT(signal,method,iterations,altermethod,nofsifts,threstype,T_mult,M1,IM2)
%EMDdenoised  : denoised signal
%signal       : Noisy signal
%method       : EMDdenoising method 
%               'IT' for Interval Thresholding (see [1]),
%               'IIT' for Iterative Interval Thresholding (see [1]),
%               'CIIT' for Clear first Iterative Interval Thresholding [1],
%iterations   : Number of averaging iterations for IIT and CIIT methods
%altermethod  : Noise altering method
%               'circ' for random circulations
%               'perm' for random permutations
%nofsifts     : Number of sifting iterations (it should take a value
%               between 5 and 10.
%threstype    : Thresholding method
%               'hard' for hard Thresholding,
%               'soft' for soft Thresholding (see [2]),
%               'softSCAD' for smoothly clipped absolute deviation
%                          (SCAD) penalty Thresholding (see [2]).
%T_mult       : Multiplication factor of the universal threshold. For
%               example, if T_mult=0.7, then the threshold applied is
%               T=0.7*sigma*sqrt(2*log*n)
%M1           : Value of parameter M1 in the reconstruction equation (see
%               [1], Eq. (11)).
%IM2          : Sets parameter M2 (see [1], Eq. (11)) equal to the number 
%               of IMFs resulted from EMD decomposition minus IM2.
% REFERECIES
% [1] Y. Kopsinis, S. McLaughlin, �Development of EMD-based Denoising
% Methods Inspired by Wavelet Thresholding,� IEEE Trans. on Signal
% Processing, VOL. 57, NO. 4, APRIL 2009.
% [2] Y. Kopsinis, S. McLaughlin, �Empirical Mode Decomposition Based
% Soft-Thresholding,� EUSIPCO 2008.

% Modified by Petteri Teikari 

    warning('off','MATLAB:dispatcher:InexactMatch')

    n=length(signal);
    t=1:n;
    if strcmp(method,'IT')==1
        iterations=1;
    end
    [IMF,localmean] = emdfull(signal,t,nofsifts);
    nofIMFs=size(IMF,1);
    clear localmean

    estimenergy_F = zeros(1,size(IMF,1)+3);
    estimenergy_F(1)=(median(abs(IMF(1,:)))/0.6745)^2;
    for k=2:size(IMF,1)+3
        estimenergy_F(k)=estimenergy_F(1)/0.719*2.01^(-k);
    end

    T_mult=T_mult*sqrt(2*log(n));

    M2=nofIMFs-IM2;
    if strcmp(method,'CIIT')==1
        clearfirst=1;
    else
        clearfirst=0;
    end

    triggered=alteringnoise(IMF,iterations,nofsifts,altermethod,clearfirst);
    EMDdenoised=EMDdenoise_averaging(IMF,triggered.IMFpros,triggered.zcpos_pr,triggered.extrema_pr,estimenergy_F,threstype,T_mult,M1,M2);

end


function out_aver=EMDdenoise_averaging(IMF,IMFpros,zcpos,extrema,IMFenergy,threstype,mult,M1,M2)
    % IMF: Matrix containing the IMFs
    % IMFpros: Cell array containing matrices of IMFs of each different random realization
    % IMFenergy: Vector containing the estimated variances of each IMF. If
    % isempty [], then the energies of the n_a processed IMFpros are estimated
    % separately. If it is equal to 0 then the average of all the IMFpros{} is
    % adopted as the final variance
    %
    % mult: The multiplication factor for the Threshold computation. T=std*mult
    % Lcount: The number up to which the first IMFs are excluded from the sum
    % Lcount_end: In the output summations, the last Lcount_end IMFs are taken
    % not from the thresholded IMFs but from the original ones
    %
    % out_aver: Cell array {Lcount,noiseaverages} containing the denoised signal that
    % corresponds to the specific i-th Lcount and averaging up to n_a, i.e.  1:n_a
    noiseaverages=length(IMFpros);

    %out = zeros(noiseaverages, ?)
    for n_a=1:noiseaverages

        out_all=EMDthres_mine(IMFpros{n_a},zcpos{n_a},extrema{n_a},IMFenergy,threstype,mult);

            out(n_a,:)=IMFsubsums(out_all,M1,IMFpros{n_a},M2);
    end
    out_aver=mean(out,1);

end

%% IMFsubsums
function out=IMFsubsums(IMFthresholded,M1,IMF,M2)
% IMFthresholded: Matrix with the thresholded IMFs.
% Lcount: The number up to which the first IMFs are excluded from the sum
% IMF: the original IMFs
% Lcount_end, the IMF number up to the end that the IMFs of the orginal
% ones are participating in the final sum.
%
% Produces a matrix, out{1,:), out(2,:), ..., out(Lcount,:), with the
% corresponding subsums of IMFthresholded. If nargin=4 then the last
% Lcount_end summed IMFs are the original and not the thresholded ones.

if nargin==2
    out=sum(IMFthresholded(M1:end,:),1);
else
    out=sum(IMFthresholded(M1:M2,:),1);
    out=out+sum(IMF(M2+1:end,:),1);
end

end
%% EMDthres_mine
function IMF_thresholded=EMDthres_mine(IMF,zcpos,extrema_Pr,IMFenergy,threstype,mult)

    % IMF: Matrix containing the IMFs
    % IMFenergy: Vector containing the estimated variances of each IMF
    % mult: The multiplication factor for the Threshold computation. T=std*mult
    extrema=extrema_Pr.extrema;
    IMFflag=extrema_Pr.IMFflag;

    nofIMFs=size(IMF,1);
    noOfSamples = size(IMF,2);
    IMF_thresholded = zeros(nofIMFs, noOfSamples);
    
    for i=1:nofIMFs
        if i<=length(IMFenergy)
            T = sqrt(IMFenergy(i))*mult;
        else
            T=0;
        end

        md=min(diff(zcpos{i}));
        if isempty(md)
            md=0;
        end
        if md<=2 || IMFflag{i}==0  %i==nofIMFs
            IMF_thresholded(i,:)=thresholdIMF_no_iterp(IMF(i,:),T,threstype,zcpos{i});
        else
            IMF_thresholded(i,:)=thresholdIMF_no_iterp2(IMF(i,:),T,threstype,zcpos{i},extrema{i});
        end
    end
end

%%
function triggered=alteringnoise(IMF,noiseaverages,maxiter,altermethod,clearfirst)

n=length(IMF);
nofIMFs=size(IMF,1);
t=1:n;

    if clearfirst==0
        noise=IMF(1,:);
        signal=sum(IMF(2:nofIMFs,:),1);
    else
        IMF1=IMF(1,:);
        qmf=MakeONFilter('Symmlet',8);
        IMF1cleaned = recdecompsh(IMF1,qmf);
        noise=IMF1-IMF1cleaned;
        signal=sum(IMF(2:nofIMFs,:),1)+IMF1cleaned;
    end

IMFflag = cell(size(IMF,1),1);
zerocrossings = zeros(1,size(IMF,1));
energyrobust = zeros(1,size(IMF,1));
for i=1:size(IMF,1)
    [maxim,minim,zc]=extremes(IMF(i,:),t);
    zcpos{i}=zc;
    extrema{i}=union(maxim,minim);
    if any(IMF(i,maxim)<0) || any(IMF(i,minim)>0)
        IMFflag{i}=0;
    else
        IMFflag{i}=1;
    end
    zerocrossings(i)=length(zc);    
    % energy(i) = cov(IMF(i,:)); % size mismatch (Petteri, not used for anything though?)    
    %energyrobust(i)=(median(abs(IMF(i,:)))/0.6745)^2; % not used either
end
zcpos_pr{1}=zcpos;
IMFpros{1}=IMF;
extrema_pr{1}.extrema=extrema;
extrema_pr{1}.IMFflag=IMFflag;
zerocrossings_pr{1}=zerocrossings;


for n_a=2:noiseaverages
    if strcmp(altermethod,'circ')==1
    noisetrig=circshift(noise',round(rand*length(noise)))';
    elseif strcmp(altermethod,'perm')==1
    noisetrig=noise(randperm(length(noise)));
    end
    
    signaltriggered=noisetrig+signal;
    [IMFtrig,localmean] = emdfull(signaltriggered,t,maxiter);
    %clear signaltriggered noisetrig
    nofIMFs_trig=size(IMFtrig,1);

    clear localmean
    
    for i=1:size(IMFtrig,1)
        [maxim,minim,zc]=extremes(IMFtrig(i,:),t);
        zcpos_tr{i}=zc;
        extrema_tr{i}=union(maxim,minim);
        if any(IMFtrig(i,maxim)<0) | any(IMFtrig(i,minim)>0)
            IMFflag_tr{i}=0;
        else
            IMFflag_tr{i}=1;
        end
        zerocrossings_tr(i)=length(zc);
    end

    IMFpros{n_a}=IMFtrig;
    clear IMFtrig
    zcpos_pr{n_a}=zcpos_tr;
    extrema_pr{n_a}.extrema=extrema_tr;
    extrema_pr{n_a}.IMFflag=IMFflag_tr;
    zerocrossings_pr{n_a}=zerocrossings_tr;
    %     end

end
triggered.IMFpros=IMFpros;
triggered.zcpos_pr=zcpos_pr;
triggered.extrema_pr=extrema_pr;

end
function [IMF,localmean] = emdfull(signal,t,maxiter);

if size(signal,2)==1, signal=signal.'; end
if size(t,2)==1, t=t'; end
nofextr = 10000;
imfnum = 1;
if imfnum==1
    axx=[min(signal),max(signal)];
end

while nofextr>2
    iter = 1;
    options.iter=iter;
    [m,siftdetails]=sifting(signal,t);
    if isempty(siftdetails)
        break
    end

    localmean(imfnum,:) = m;

    h = signal - m;
    if maxiter>=iter
        crettotal = 1;
    else
        crettotal = 0;
    end

    iter = iter+1;

    while crettotal
        [m,siftdetails]=sifting(h,t);
        if isempty(siftdetails)
            break
        end

        if maxiter>=iter
            crettotal = 1;
        else
            crettotal = 0;
        end

        if crettotal==1
            h = h - m;
            localmean(imfnum,:) = localmean(imfnum,:) + m;
        end

        iter = iter+1;
    end
    IMF(imfnum,:) = h;
    signal = signal - h;
    imfnum = imfnum + 1;
end
IMF(imfnum,:)=signal;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m,siftdetails]=sifting(signal,t)
interporder=3;
[maxi, mini, zerocross] = extremes(signal,t);

nofextr = length(maxi)+length(mini);

if length(maxi) >= 2 && length(mini) >= 2
    [maxp_ex, maxv_ex, minp_ex, minv_ex] = edge_extrapolation(maxi, mini, signal, t, 3, 'mirror');
    if length(maxp_ex) <= interporder || length(minp_ex) <= interporder
        m=zeros(1,length(signal));
        siftdetails=[];
        return
    end
else
    m=zeros(1,length(signal));
    siftdetails=[];
    return
end


envmax = interp1(maxp_ex,maxv_ex,t,'spline'); % envelop of the max
envmin = interp1(minp_ex,minv_ex,t,'spline'); % envelop of the min



m = 1/2*(envmax + envmin);
siftdetails.envmax=envmax;
siftdetails.envmin=envmin;
siftdetails.maxi=maxi;
siftdetails.mini=mini;
siftdetails.zerocross=zerocross;
siftdetails.nofextr=nofextr;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function [maxima , minima, zerocross] = extremes(signal,t)

dif = diff(signal); %approximate derivative
%   h=t(2)-t(1);
%   dif=nd5p(signal,h,length(signal));
dif = sign(dif);

% It considers the appearence of single zeros only.
if any(dif==0)
    zdif=find(dif==0);      %in case of zero the sign of the previous
    if zdif(1)==1, zdif(1)=[]; end
    if zdif(end)==length(dif), zdif(end)=[]; end
    dif(zdif)=dif(zdif-1);  %point is adopted.
end

N = length(dif);
%sdif = dif(1:N-1)+ dif(2:N);
%extremespos = find(sdif==0);

sdif = dif(1:N-1).*dif(2:N);
extremespos = find(sdif==-1);
extremestype = dif(extremespos);

maxima = extremespos(extremestype>0)+1;
minima = extremespos(extremestype<0)+1;

zerocross = find(signal(1:N-1).*signal(2:N)<0);

end

function [maxp_ex, maxv_ex, minp_ex, minv_ex] = edge_extrapolation(maxi, mini, signal, t, N_extr, method);

% This function contains some parts (possibly modified) from the EMD.m function made by
%
% G. Rilling, July 2002
%
% that computes EMD (Empirical Mode Decomposition) according to:
%
% N. E. Huang et al., "The empirical mode decomposition and the
% Hilbert spectrum for non-linear and non stationary time series analysis,"
% Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998
%
% with variations reported in:
%
% G. Rilling, P. Flandrin and P. Gon?alv?s
% "On Empirical Mode Decomposition and its algorithms"
% IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
% NSIP-03, Grado (I), June 2003




% Extrapolates a number of N_extr extemes around the edges of the signal
%
% method: 'mirror' referes to the method of Rilling et al., "on empirical
%          mode decomposition and its algorithms",
%          sortly, new extrapolation methods will be added.
%
% maxi (mini): the index positions (1,2,...,N) of the maxima (minima)
% N_extr: The number of extrapolating points

% maxp_ex (minp_ex): the time positions of the maxima (minima) including
%                    the extrapolated points at the edges
% maxv_ex (minv_ex): the values of the maxima (minima) including
%                    the values of the extrapolated points at the edges

maxp = t(maxi); % maxp (minp): the time positions of the maxima (minima)
minp = t(mini);
maxv = signal(maxi); % maxv (minv): the values of the maximuma (minimuma)
minv = signal(mini);


if strcmp(method,'mirror')==1
    indmax = maxi;
    indmin = mini;
    m = signal;
    NBSYM = N_extr;
    lx = length(signal);
    %% copy from the Rilling - Flandrin m-file.
    if indmax(1) < indmin(1)
        if m(1) > m(indmin(1))
            lmax = fliplr(indmax(2:min(end,NBSYM+1)));
            lmin = fliplr(indmin(1:min(end,NBSYM)));
            lsym = indmax(1);
        else
            lmax = fliplr(indmax(1:min(end,NBSYM)));
            lmin = [fliplr(indmin(1:min(end,NBSYM-1))),1];
            lsym = 1;
        end
    else
        if m(1) < m(indmax(1))
            lmax = fliplr(indmax(1:min(end,NBSYM)));
            lmin = fliplr(indmin(2:min(end,NBSYM+1)));
            lsym = indmin(1);
        else
            lmax = [fliplr(indmax(1:min(end,NBSYM-1))),1];
            lmin = fliplr(indmin(1:min(end,NBSYM)));
            lsym = 1;
        end
    end

    if indmax(end) < indmin(end)
        if m(end) < m(indmax(end))
            rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
            rmin = fliplr(indmin(max(end-NBSYM,1):end-1));
            rsym = indmin(end);
        else
            rmax = [lx,fliplr(indmax(max(end-NBSYM+2,1):end))];
            rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
            rsym = lx;
        end
    else
        if m(end) > m(indmin(end))
            rmax = fliplr(indmax(max(end-NBSYM,1):end-1));
            rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
            rsym = indmax(end);
        else
            rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
            rmin = [lx,fliplr(indmin(max(end-NBSYM+2,1):end))];
            rsym = lx;
        end
    end

    tlmin = 2*t(lsym)-t(lmin);
    tlmax = 2*t(lsym)-t(lmax);
    trmin = 2*t(rsym)-t(rmin);
    trmax = 2*t(rsym)-t(rmax);

    % in case symmetrized parts do not extend enough
    if tlmin(1) > t(1) | tlmax(1) > t(1)
        if lsym == indmax(1)
            lmax = fliplr(indmax(1:min(end,NBSYM)));
        else
            lmin = fliplr(indmin(1:min(end,NBSYM)));
        end
        if lsym == 1
            error('bug')
        end
        lsym = 1;
        tlmin = 2*t(lsym)-t(lmin);
        tlmax = 2*t(lsym)-t(lmax);
    end

    if trmin(end) < t(lx) | trmax(end) < t(lx)
        if rsym == indmax(end)
            rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
        else
            rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
        end
        if rsym == lx
            error('bug')
        end
        rsym = lx;
        trmin = 2*t(rsym)-t(rmin);
        trmax = 2*t(rsym)-t(rmax);
    end

    mlmax =m(lmax);
    mlmin =m(lmin);
    mrmax =m(rmax);
    mrmin =m(rmin);


    if length(mlmax)<NBSYM
        KK=NBSYM-length(mlmax);
        if length(mlmax)>1
            dml=mlmax(1)-mlmax(2);
            dtl=abs(tlmax(1)-tlmax(2));
        else
            dml=mlmax-maxv(1);
            dtl=abs(tlmax-maxp(1));
        end
        mlmax=[fliplr(mlmax(1)+dml*(1:KK)) mlmax];
        tlmax=[fliplr(tlmax(1)-dtl*(1:KK)) tlmax];
    end


    if length(mrmax)<NBSYM
        KK=NBSYM-length(mrmax);
        if length(mrmax)>1
            dmr=mrmax(end)-mrmax(end-1);
            dtr=abs(trmax(end)-trmax(end-1));
        else
            dmr=mrmax-maxv(end);
            dtr=abs(trmax-maxp(end));
        end

        mrmax=[mrmax mrmax(end)+dmr*(1:KK)];
        trmax=[trmax trmax(end)+dtr*(1:KK)];
    end


    if length(mlmin)<NBSYM
        KK=NBSYM-length(mlmin);
        if length(mlmin)>1
            dml=mlmin(1)-mlmin(2);
            dtl=abs(tlmin(1)-tlmin(2));
        else
            dml=mlmin-minv(1);
            dtl=abs(tlmin-minp(1));
        end
        mlmin=[fliplr(mlmin(1)+dml*(1:KK)) mlmin];
        tlmin=[fliplr(tlmin(1)-dtl*(1:KK)) tlmin];
    end


    if length(mrmin)<NBSYM
        KK=NBSYM-length(mrmin);
        if length(mrmin)>1
            dmr=mrmin(end)-mrmin(end-1);
            dtr=abs(trmin(end)-trmin(end-1));
        else
            dmr=mrmin-minv(end);
            dtr=abs(trmin-minp(end));
        end
        mrmin=[mrmin mrmin(end)+dmr*(1:KK)];
        trmin=[trmin trmin(end)+dtr*(1:KK)];
    end
    maxv_ex = [mlmax maxv mrmax];
    maxp_ex = [tlmax maxp trmax];
    minv_ex = [mlmin minv mrmin];
    minp_ex = [tlmin minp trmin];

else
    fprintf('Unknown extrapolation method')
end

end


%% thresholdIMF_no_iterp
function thresholded=thresholdIMF_no_iterp(IMF,T,threstype,zcpos)
% It does not uses the final interpolation.

n=length(IMF);
temp_out=zeros(1,n);

temp_out(IMF > T)=1;
DD=diff([0 temp_out]);
signalpos_plus=find(DD==1);
temp_out=zeros(1,n);
temp_out(IMF < -T)=1;
DD=diff([0 temp_out]);
signalpos_minus=find(DD==1);
signalpos=sort([signalpos_plus signalpos_minus]);

temp_out=zeros(1,n);
zerocrossings=[0 zcpos n]+0.5;
nodesexclude=[];
for fsp=1:length(signalpos)
    ff=find(zerocrossings<signalpos(fsp));
    ff=length(ff);
    if zerocrossings(ff)~=n
        if strcmp(threstype,'hard')==1
            temp_out(ceil(zerocrossings(ff)):floor(zerocrossings(ff+1)))=1;
        elseif strcmp(threstype,'soft')==1
            pick=max(abs(IMF((ceil(zerocrossings(ff)):floor(zerocrossings(ff+1))))));
            temp_out(ceil(zerocrossings(ff)):floor(zerocrossings(ff+1)))=(pick-T)/pick;
        elseif strcmp(threstype,'softSCAD')==1
            pick=max(abs(IMF((ceil(zerocrossings(ff)):floor(zerocrossings(ff+1))))));
            alpha=3.7;
            if pick<=2*T
                temp_out(ceil(zerocrossings(ff)):floor(zerocrossings(ff+1)))=max(0,(pick-T))/pick;
            elseif pick<=alpha*T
                temp_out(ceil(zerocrossings(ff)):floor(zerocrossings(ff+1)))=(((alpha-1)*pick-alpha*T)/(alpha-2))/pick;
            else
                temp_out(ceil(zerocrossings(ff)):floor(zerocrossings(ff+1)))=1;
            end
        end
    end
end
temp_out1=IMF.*temp_out;
thresholded=temp_out1;
end

%% thresholdIMF_no_iterp2
function thresholded=thresholdIMF_no_iterp2(IMF,T,threstype,zcpos,extremapos)
% It does not uses the final interpolation.

if length(extremapos)>length(zcpos)+1 | length(zcpos)<=1 %The trend IMF
    thresholded=IMF;
else

    n=length(IMF);

    not_thresholded=find(abs(IMF(extremapos)) >= T);


    if extremapos(1)<=zcpos(1)
        if extremapos(2)<=zcpos(1)
            zcpos=[1 round((extremapos(2)+extremapos(1))/2) zcpos];
        else
            zcpos=[1 zcpos];
        end
    end
    if extremapos(end)>=zcpos(end)
        if extremapos(end-1)>=zcpos(end)
            zcpos=[zcpos round((extremapos(end)+extremapos(end-1))/2) n];
        else
            zcpos=[zcpos n];
        end
    end

    temp_out=zeros(1,n);
    for fsp=1:length(not_thresholded)
        if strcmp(threstype,'hard')==1
            temp_out(zcpos(not_thresholded(fsp))+1:zcpos(not_thresholded(fsp)+1))=1;
        elseif strcmp(threstype,'soft')==1
            pick=max(abs(IMF(zcpos(not_thresholded(fsp))+1:zcpos(not_thresholded(fsp)+1))));
            temp_out(zcpos(not_thresholded(fsp))+1:zcpos(not_thresholded(fsp)+1))=(pick-T)/pick;

        elseif strcmp(threstype,'softSCAD')==1
            pick=max(abs(IMF(zcpos(not_thresholded(fsp))+1:zcpos(not_thresholded(fsp)+1))));
            alpha=3.7;
            if pick<=2*T
                temp_out(zcpos(not_thresholded(fsp))+1:zcpos(not_thresholded(fsp)+1))=max(0,(pick-T))/pick;
            elseif pick<=alpha*T
                temp_out(zcpos(not_thresholded(fsp))+1:zcpos(not_thresholded(fsp)+1))=(((alpha-1)*pick-alpha*T)/(alpha-2))/pick;
            else
                temp_out(zcpos(not_thresholded(fsp))+1:zcpos(not_thresholded(fsp)+1))=1;
            end
        end
    end

    % In order to let the fist and last parts of IMF if they are over the
    % threshold
    if abs(IMF(1)) >= T;
        temp_out(1:zcpos(1))=1;
    end
    if abs(IMF(end)) >= T;
        temp_out(zcpos(end)+1:n)=1;
    end

    thresholded=IMF.*temp_out;
end
end