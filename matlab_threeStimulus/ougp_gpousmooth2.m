% From: http://clinengnhs.liv.ac.uk/
function [z,nf]=ougp_gpousmooth2(x,y,wc,type)

    % Analytical form of inverse covariance for a generalised
    % 1-dimensional Ornstein-Uhlenbeck process

    %
    % INPUTS
    % x - times
    % y - data
    % wc - cutoff frequency (or 2-vector of frequencies for bandpass)
    %      (expressed as fraction of Nyquist frequency)  

    % type - string: lp, hp, bp (for low-pass, high-pass or band-pass)
    %       default: lp

    % OUTPUTS
    % z - smoothed data
    % nf - estimate of the Nyquist frequency

    % A Eleuteri - 27/02/2009
    % -----------------------------------------------------------------------
    % Reference paper: G. B. Rybicki and W. H. Press, Physical Review Letters
    % Vol. 74, Num. 7, p.1060-1063
    % http://dx.doi.org/10.1007/s11517-012-0928-2

    if nargin<4
        type='lp';
    end

    if length(wc)==2
        type='bp';
    end

    nf=median(1./diff(x))/2; % estimate of the Nyquist frequency

    if (strcmp(type,'lp'))
        z=gouinvcovfilt(x,y,wc*nf,type);
    elseif strcmp(type,'hp')
        z=gouinvcovfilt(x,y,wc*nf,type);
    elseif strcmp(type,'bp')
        tmp=gouinvcovfilt(x,y,wc(2)*nf,'lp');
        z=gouinvcovfilt(x,tmp,wc(1)*nf,'hp');
    else
        error('unrecognised option or wrong number of frequencies.');
    end

    return



    function z=gouinvcovfilt(x,y,wc,type)
        if length(wc)>1
            error(['frequency must be a scalar for the ',type,' option.']);
        end

        T=length(x);

        if (strcmp(type,'lp'))
            k=sqrt(2)*pi*(sqrt(2)-1)^(-1/4);
        else
            k=sqrt(2)*pi*(sqrt(2)-1)^(1/4);
        end

        w=k*(1+1i)*wc*diff(x);
        r=exp(-w);
        e=r./(1-r.^2);
        re=r.*e;
        d=1+[re;0]+[0;re];

        invK=spdiags([-[e;0],d,-[0;e]],-1:1,T,T); % inverse covariance matrix
        f=cmplxtrf(y,w);
        u=(invK)\(f);

        if (strcmp(type,'lp'))
            z=y-real(u);
        else
            z=real(u);
        end
        
        return



    function f=cmplxtrf(y,w)

        f=zeros(size(y));
        T=length(y);
        dy=diff(y);
        f(1)=-0.5*dy(1)/w(1);
        f(2:T-1)=0.5*(-dy(2:end)./w(2:end)+dy(1:end-1)./w(1:end-1));
        f(T)=0.5*dy(end)/w(end);

    return
