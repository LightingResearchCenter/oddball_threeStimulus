function Example

SNRdB=5;
n=2048;

x = MakeSignal('Doppler',n);
t=1:length(x);

%t = linspace(-pi,pi,n);
%x = cos(t)

sigpower=mean(x.^2);
noisepower=sigpower/(10^(SNRdB/10));
noisesignal=randn(1,length(x))*sqrt(noisepower);
signal=x+noisesignal;

% plot(t,signal,'b', t,noisesignal, 'r', t, x, 'k', 'LineWidth', 2)

n_a=5;
M1=3;
IM2=2;
T_mult=0.7;
threstype='soft';
nofsifts=8;
iterations=15;
method='CIIT';
altermethod='perm';
[EMDdenoised, IMF] = EMDdenoise_modPT(signal,method,iterations,altermethod,nofsifts,threstype,T_mult,M1,IM2);

meanx2=mean(x.^2);
SNR=10*log10(meanx2/mean((x-EMDdenoised).^2));

qmf=MakeONFilter('Symmlet',8);
out_A1 = recdecompsh(signal,qmf);
SNR_A1=10*log10(meanx2/mean((x-out_A1).^2));   

figure;
plot(signal, 'Color', [0 .55 1])
hold on
plot(x, 'k')
plot(EMDdenoised,'r--')
plot(out_A1,'g--')

    legend('Input with Noise', 'Input IDEAL', 'EMDdenoised', 'Reconstruction?','Location', 'Best')
        legend('boxoff')
    xlim([min(t) max(t)])
    


end


