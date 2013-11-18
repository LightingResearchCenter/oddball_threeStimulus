function [spec_wt,scale] = plot_wavespec_cwt(signal,sc_numb)

% Plot function of the scale continuous wavelet spectrum
%
% spec_wt = plot_wavespec_cwt(signal,sc_numb)
%
% INPUT PARAMETER---------------------------------------
%
% signal:       Input data series
% sc_numb:      Number of scales in time-scale plot;
%
% OUTPUT VARIABLES--------------------------------------
% 
% spec_wt:      [sc_numb,length(signal)] matrix of the time dependent spectral power
%
% EXAMPLE-----------------------------------------------
%
% N=1024;
% t=1:N;
% Ht=0.5+0.5.*(sin(0.0025.*pi.*t));
% [mBm,mGn]=mBm_mGn(N,Ht);
% [spec_wt, scale] = plot_wavespec_cwt(mGn,100);
%
% ------------------------------------------------------
% written be Espen A. F. Ihlen, 2010 (espen.ihlen@ntnu.no)

N=length(signal);
t=1:N;
scale=linspace(1/6.5,(N/2)/6.5,sc_numb);
wt = cwt(signal,scale,'morl');
scale=linspace(1,N/2,sc_numb);
spec_wt=(abs(wt)).^2;
figure;
subplot(211)
plot(t,signal);
xlabel('t');ylabel('amplitude');title('time series')
axis([1 1024 -5 5])
subplot(212)
mesh(t,scale,spec_wt);
view([0 90]);axis([1 N 0 500])
%colormap(1-gray(256))
xlabel('t');ylabel('\Deltat');title('Wavelet spectrum')