% Modified from HRVAS toolbox
function [cwtpower,WT,f,scale,Cdelta,n,dj,dt,variance,coi]= analyze_waveletWrapper(y,fs,timeDiv,k)
    
    variance = std(y)^2;
    y = (y - mean(y)) / sqrt(variance);

    n = length(y);
    dt = 1/fs;   
    %xlim = [0,t2(end)];  % plotting range
    pad = 1;      % pad the time series with zeroes (recommended)
    dj = 1/64;    % this will do 4 sub-octaves per octave
    s0 = timeDiv*dt;    % e.g. 2*dt (half the sampling frequency)
    j1 = 7/dj;    % this says do 7 powers-of-two with dj sub-octaves each
    lag1 = 0.72;  % lag-1 autocorrelation for red noise background

    mother = 'Morlet';
    %mother = 'DOG';
    %mother = 'Paul';
    Cdelta = 0.776;   % this is for the MORLET wavelet
    
    % Wavelet transform
    [wave,period,scale,coi] = analyze_waveletTransform(y,dt,pad,dj,s0,j1,mother,k);
        % Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
        % Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.    
        
        % coi = if specified, then return the Cone-of-Influence, which is a vector
        %        of N points that contains the maximum period of useful information
        %        at that particular time.
        %        Periods greater than this are subject to edge effects.
        %        This can be used to plot COI lines on a contour plot by doing:
        %
        %              contour(time,log(period),log(power))
        %              plot(time,log(coi),'k')
    
    cwtpower = (abs(wave)).^2 ; % compute wavelet power spectrum              
    %cwtpower = wave; % compute wavelet power spectrum              
    
    % match the variable name to ERPWAVELAB
    WT = wave;
    
    f=fliplr(1./period); %frequency in ascending order
    cwtpower=flipud(cwtpower); %flip to match freq. order
 