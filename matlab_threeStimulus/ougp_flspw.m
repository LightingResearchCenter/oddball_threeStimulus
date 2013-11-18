% From: http://clinengnhs.liv.ac.uk/
function [Pr,fr,A,z0,A0,ofac]=ougp_flspw(x,y,f,p0,iofac)

    % The Fastest Lomb-Scargle Periodogram estimator in the West
    % [Pr,fr,A,z0,A0,ofac]=fLSPw(x,y,f,p0,iofac)

    %  This is a very efficient MATLAB implementation of the Press-Rybicki ...
    %     approximation of the Lomb-Scargle's periodogram estimator.

    %
    % INPUT
    % x -   time vector
    % y -   amplitude vector
    % f -   frequency vector which define the frequency bands of interest. Try to
    %           avoid very small lower bounds (i.e. <1e-4) to reduce the computational
    %           complexity of the algorithm.
    % p0 -  false alarm probability for detection (default=0.05)
    % iofac - oversampling factor (default=4)
    %

    % OUTPUT
    % Pr -  cell array of periodogram estimates for each frequency band.
    %           NOTE: divide the periodogram estimate by the sampling frequency to have a
    %           scale compatible with MATLAB's pwelch routine.
    % fr -  cell array of frequencies for each frequency band
    % A -   vector of integrals of the psd for each frequency band
    % z0 -  vector of thresholds on psd amplitude for detection of signals in frequency band
    % A0 -  vector of thresholds on integrals of psd for detection of signals in
    %       frequency band
    % ofac -    effective oversampling factor, different from input when the
    %           lowest frequency required is smaller than the one defined by iofac

    % A Eleuteri 23/7/2008, 16/12/2008

    if nargin<5
        iofac=4;
    end

    if nargin<4
        p0=0.05;
    end

    flowbound=1e-4; % lower bound on frequencies, decrease only if you know that ...
                    % you want to analyse data at very small frequencies

    % constants (usually should not need change them)
    nfreq=64;
    macc=4; % accuracy of approximation

    f=sort(f); % list of frequencies which define the frequency bands of interest
    fmin=max(f(1),flowbound); % minimum frequency
    n=length(y);
    ave=mean(y);
    vr=var(y);
    xmin=min(x);
    xmax=max(x);
    xdif=xmax-xmin;
    ofac=max(iofac,1/(xdif*fmin));  % if the lowest frequency is not compatible with ...
                                    % the default oversampling factor, redefine it
    df=1/(xdif*ofac);
    fc=n/(2*xdif); % "average" Nyquist frequency

    hifac=f/fc;
    nout=fix(0.5*ofac*hifac*n);
    noutmax=nout(end);
    nfreqt=2*noutmax*macc;
    %nfreq=nfreqt;

    if nfreq<nfreqt
        nfreq=2^nextpow2(nfreqt);
        %nfreq=2*nfreq;
    end

    %ndim=2*nfreq;
    ndim=nfreq;

    % "extirpolate" the data
    wk1=zeros(ndim,1);
    wk2=wk1;
    fac=ndim/(xdif*ofac);
    fndim=ndim;
    ck=1+rem((x-xmin)*fac,fndim);
    ckk=1+rem(2*(ck-1),fndim);

    for j=1:n
       wk1=spread(wk1,y(j)-ave,ndim,ck(j),macc);
       wk2=spread(wk2,1,ndim,ckk(j),macc);
    end

    tmp1=fft(wk1(1:nfreq)');
    tmp2=fft(wk2(1:nfreq)');

    % Put the fft output in a format compatible with the NR implementation
    wk1=fliplr(tmp1(length(tmp1)/2+2:end));
    wk1=[real(wk1);imag(wk1)];
    wk1=wk1(:);
    wk2=fliplr(tmp2(length(tmp2)/2+2:end));
    wk2=[real(wk2);imag(wk2)];
    wk2=wk2(:);

    % Compute the Lomb-Scargle periodogram
    k=1:2:(2*noutmax-1);
    kp1=k+1;
    hypo=sqrt(wk2(k).^2+wk2(kp1).^2);
    hc2wt=0.5*wk2(k)./hypo;
    hs2wt=0.5*wk2(kp1)./hypo;
    cwt=sqrt(0.5+hc2wt);
    swt=abs(sqrt(0.5-hc2wt)).*sign(hs2wt);
    den=0.5*n+hc2wt.*wk2(k)+hs2wt.*wk2(kp1);
    cterm=(cwt.*wk1(k)+swt.*wk1(kp1)).^2./den;
    sterm=(cwt.*wk1(kp1)-swt.*wk1(k)).^2./(n-den);
    P=(cterm+sterm)/(2*vr);
    F=(1:noutmax)'*df;

    % Store the periodogram segments corresponding to the frequency ranges
    fr=cell(length(f)-1,1);
    Pr=fr;

    for n=1:(length(f)-1)

        fr{n}=F((nout(n)+1):nout(n+1));
        Pr{n}=P((nout(n)+1):nout(n+1));

        % evaluate the total power in the band (f0,f1] by rectangular integration
        A(n)=df*sum(Pr{n});  

        % detection statistics
        z0(n)=log(length(fr{n}))-log(p0);   % approximate threshold for spectral density peaks. ...
                                            % Assumes p0<<1

        A0(n)=gaminv(1-p0,length(fr{n}),df); % threshold for integrated power over band

    end

% Extirpolation (See "Numerical Recipes" by Press et al.)
function yy=spread(yy,y,n,x,m)

        nfac=[1,1,2,6,24,120,720,5040,40320,362880];
        ix=x;

        if(x==single(x))
            yy(ix)=yy(ix)+y;
        else
        ilo=min(max( round(x-0.5*m+1),1),  n-m+1);
        ihi=ilo+m-1;
        nden=nfac(m);
        fac=x-ilo;

        j=(ilo+1):ihi;
        fac=fac*cumprod(x-j);fac=fac(end);
        yy(ihi)=yy(ihi)+y*fac/(nden*(x-ihi));   

        for j=(ihi-1):-1:ilo
            nden=fix(nden/(j+1-ilo))*(j-ihi);
            yy(j)=yy(j)+y*fac./(nden*(x-j));
        end

    end



%%%%%%%% end of code
