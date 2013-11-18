% Generate auditory stimuli for auditory oddball task
function generateAudioWav_threeStimulus()

    % Petteri Teikari, 2013
    % petteri.teikari@gmail.com
    % Lighting Research Center, Troy, NY, USA    
    close all
    scrsz = get(0,'ScreenSize'); % get screen size for plotting

    %% SETTINGS
    
        % Styling
        style.fontName = 'latin modern roman';
        style.fontBaseSize = 9;        
       
        % FIGURE autosave options
        style.imgOutRes       = '-r300'; % dpi
        style.imgOutAntiAlias = '-a1';   % '-a0' least, '-a4' maximum
        style.imgOutautoSavePlot = 1;
    
    %% AUDIO PARAMETERS
    
        % Sampling rate
        Fs = 48000; % [Hz]
        bits = 16; % bit depth
    
        % Sinusoidal frequencies, see e.g. Frank et al. (2012)
        % http://dx.doi.org/10.1016/j.ijpsycho.2012.04.005
        standard.f = 500; % [Hz]
        oddball.f  = 1000; % [Hz] (target)
        
        % White Noise (non-target / distracter)
        % bandwidth-limited "hissing sound", so by defnition not really
        % white anymore but keeping the same definitions as used in Frank et al. (2012)
        noise.f = [500 5000]; % [lowCut hiCut]    
    
        % Audio durations
        standard.t = 0.060; % [s]
        oddball.t = standard.t;
        noise.t = standard.t;
        
        % Audio amplitudes
        standard.volume = 0.8; % normalized to something, will clip at 1.0
        oddball.volume  = standard.volume;
        noise.volume    = standard.volume;
        
        % SOA/ISI between audio stimuli,
        % Strictly speaking the time between the audio offset, 
        % and the onset of next audio stimulus
        SOA = 1.930; % [s] + standard.t
        
        % Rise and fall time
        riseTime = 0.005; % [s]
        fallTime = 0.005; % [s]
    
    
    %% GENERATE    
    
        standard.noOfSamples = standard.t * Fs; % [s * samples/s = samples] % e.g. 9600 samples
        oddball.noOfSamples = oddball.t * Fs; % [s * samples/s = samples] % e.g. 9600 samples
        noise.noOfSamples = noise.t * Fs;
        
        standard.timeVector = (linspace(0,standard.t, standard.noOfSamples))';
        oddball.timeVector = (linspace(0,oddball.t, oddball.noOfSamples))';
        noise.timeVector = (linspace(0,noise.t, noise.noOfSamples))';
        
        standard.noOfCycles = standard.noOfSamples / standard.f;
        oddball.noOfCycles  = oddball.noOfSamples / oddball.f;           
        
        standard.wave = standard.volume * sin(2*pi*standard.timeVector * standard.f);
        oddball.wave = oddball.volume * sin(2*pi*oddball.timeVector * oddball.f);
        
        % Generate the white noise
        noise.wave = oddball.volume * generateWhiteNoiseWithBandwidthLimitation(noise.timeVector, noise.f, Fs);
        
            % save for calibration purposes without the envelope and the
            % silent inter-stimulus 
            fileName = 'standardTone_calibration.wav';
            try
                audiowrite(fileName, standard.wave, Fs, 'BitsPerSample', bits)
            catch err
                err
                wavwrite(standard.wave, Fs, bits, fileName)
            end
            
            fileName = 'targetTone_calibration.wav';
            try
                audiowrite(fileName, oddball.wave, Fs, 'BitsPerSample', bits)
            catch err
                err
                wavwrite(oddball.wave, Fs, bits, fileName)
            end
            
            fileName = 'noiseDistracter_calibration.wav';
            try
                audiowrite(fileName, noise.wave, Fs, 'BitsPerSample', bits)
            catch err
                err
                wavwrite(noise.wave, Fs, bits, fileName)
            end
            
    
    %% Add the rise/fall time envelopes 
    
        % See for e.g. http://dx.doi.org/10.1016/j.brainres.2008.11.087 (Fig. 7 or Fig. 8)
        % So the envelope need to be just linear ramps
        noOfSamplesPerRise = riseTime * Fs;
        noOfSamplesPerFall = fallTime * Fs;
        
        riseRamp = (linspace(0, standard.volume, noOfSamplesPerRise))';
        fallRamp = flipud((linspace(0, standard.volume, noOfSamplesPerFall))');
        
        envelope = ones(length(standard.wave),1);
            envelope(1:length(riseRamp)) = riseRamp;
            envelope(end-(length(fallRamp)-1):end) = fallRamp;
        
        standard.waveWithEnvelope = standard.wave .* envelope;
        oddball.waveWithEnvelope  = oddball.wave  .* envelope;
        noise.waveWithEnvelope  = noise.wave  .* envelope;
       
        
    %% ADD the SOA/ISI
    
        standard.tFinal = standard.t + SOA;
        oddball.tFinal  = oddball.t + SOA;
        noise.tFinal  = noise.t + SOA;
    
        standard.noOfSamplesFinal = standard.tFinal * Fs;
        oddball.noOfSamplesFinal  = oddball.tFinal * Fs;
        noise.noOfSamplesFinal  = noise.tFinal * Fs;
    
        standard.finalTimeVector =(linspace(0, standard.tFinal, standard.noOfSamplesFinal))';
        oddball.finalTimeVector = (linspace(0, oddball.tFinal, oddball.noOfSamplesFinal))';
        noise.finalTimeVector = (linspace(0, noise.tFinal, noise.noOfSamplesFinal))';
        
        standard.waveWithSOA = zeros(length(standard.finalTimeVector),1);
            standard.waveWithSOA(1:length(standard.waveWithEnvelope)) = standard.waveWithEnvelope;
            
        oddball.waveWithSOA  = zeros(length(oddball.finalTimeVector),1);
            oddball.waveWithSOA(1:length(oddball.waveWithEnvelope)) = oddball.waveWithEnvelope;
            
        noise.waveWithSOA  = zeros(length(noise.finalTimeVector),1);
            noise.waveWithSOA(1:length(noise.waveWithEnvelope)) = noise.waveWithEnvelope;

    
    %% DO FFT for debugging
    
        % Without SOA, and without envelope
        
            [~, standard.waveFFT] = computeSingleSidedAmplitudeSpectrum(standard.timeVector, standard.wave, Fs);
            [~, oddball.waveFFT] = computeSingleSidedAmplitudeSpectrum(oddball.timeVector, oddball.wave, Fs);
            [~, noise.waveFFT] = computeSingleSidedAmplitudeSpectrum(noise.timeVector, noise.wave, Fs);       

        % Without SOA, and with envelope
        
            [~, standard.waveFFTWithEnvelope] = computeSingleSidedAmplitudeSpectrum(standard.timeVector, standard.waveWithEnvelope, Fs);
            [~, oddball.waveFFTWithEnvelope] = computeSingleSidedAmplitudeSpectrum(oddball.timeVector, oddball.waveWithEnvelope, Fs);
            [freqVector, noise.waveFFTWithEnvelope] = computeSingleSidedAmplitudeSpectrum(noise.timeVector, noise.waveWithEnvelope, Fs);             
                
        % With SOA
            
            [~, standard.waveFFTWithSOA] = computeSingleSidedAmplitudeSpectrum(standard.timeVector, standard.waveWithSOA, Fs);
            [~, oddball.waveFFTWithSOA] = computeSingleSidedAmplitudeSpectrum(oddball.timeVector, oddball.waveWithSOA, Fs);
            [freqVectorWithSOA, noise.waveFFTWithSOA] = computeSingleSidedAmplitudeSpectrum(noise.timeVector, noise.waveWithSOA, Fs);            
        
    
    %% PLOT
    
        fig = figure('Color', 'w');
            set(fig, 'Position', [0.5*scrsz(3) 0.1*scrsz(4) 0.48*scrsz(3) 0.8*scrsz(4)])
            
            % subplot layout
            rows = 6;
            cols = 3;

            i = 1;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(standard.timeVector, standard.wave, 'r');
                tit(i) = title(['STANDARD, f = ', num2str(standard.f), ' Hz']);

            i = 2;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(oddball.timeVector, oddball.wave, 'b');
                tit(i) = title(['ODDBALL, f = ', num2str(oddball.f), ' Hz']);
                
            i = 3;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(noise.timeVector, noise.wave, 'b');
                tit(i) = title(['NOISE, f = ', num2str(noise.f), ' Hz']);

            i = 4;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(standard.timeVector, standard.waveWithEnvelope, 'r');
                tit(i) = title(['STANDARD w Envelope, f = ', num2str(standard.f), ' Hz']);

            i = 5;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(oddball.timeVector, oddball.waveWithEnvelope, 'b');
                tit(i) = title(['ODDBALL w Envelope, f = ', num2str(oddball.f), ' Hz']);
                
            i = 6;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(noise.timeVector, noise.waveWithEnvelope, 'b');
                tit(i) = title(['NOISE w Envelope, f = ', num2str(noise.f), ' Hz']);

            i = 7;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(standard.finalTimeVector, standard.waveWithSOA, 'r');
                tit(i) = title(['STANDARD w SOA, f = ', num2str(standard.f), ' Hz']);

            i = 8;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(oddball.finalTimeVector, oddball.waveWithSOA, 'b');
                tit(i) = title(['ODDBALL w SOA, f = ', num2str(oddball.f), ' Hz']);
                
            i = 9;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(noise.finalTimeVector, noise.waveWithSOA, 'b');
                tit(i) = title(['NOISE w SOA, f = ', num2str(noise.f), ' Hz']);

            i = 10;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(freqVector, standard.waveFFT, 'r');
                [v,ind] = max(standard.waveFFT);
                freqPeak = freqVector(ind);
                tit(i) = title(['STANDARD, f = ', num2str(freqPeak,5), ' Hz']);
                    xlim([0 1.2*(max([oddball.f standard.f noise.f]))])

            i = 11;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(freqVector, oddball.waveFFT, 'b');
                [v,ind] = max(oddball.waveFFT);
                freqPeak = freqVector(ind);
                tit(i) = title(['ODDBALL, f = ', num2str(freqPeak,5), ' Hz']);
                    xlim([0 1.2*(max([oddball.f standard.f noise.f]))])
                    
            i = 12;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(freqVector, noise.waveFFT, 'b');
                [v,ind] = max(noise.waveFFT);
                freqPeak = freqVector(ind);
                tit(i) = title(['NOISE, f = ', num2str(freqPeak,5), ' Hz']);
                    xlim([0 1.2*(max([oddball.f standard.f noise.f]))])

            i = 13;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(freqVector, standard.waveFFTWithEnvelope, 'r');
                [v,ind] = max(standard.waveFFTWithEnvelope);
                freqPeak = freqVector(ind);
                tit(i) = title(['STANDARD w Envelope, f = ', num2str(freqPeak,5), ' Hz']);
                    xlim([0 1.2*(max([oddball.f standard.f noise.f]))])            
            
            i = 14;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(freqVector, oddball.waveFFTWithEnvelope, 'b');
                [v,ind] = max(oddball.waveFFTWithEnvelope);
                freqPeak = freqVector(ind);
                tit(i) = title(['ODDBALL w Envelope, f = ', num2str(freqPeak,5), ' Hz']);
                    xlim([0 1.2*(max([oddball.f standard.f noise.f]))])
            
            i = 15;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(freqVector, noise.waveFFTWithEnvelope, 'b');
                [v,ind] = max(noise.waveFFTWithEnvelope);
                freqPeak = freqVector(ind);
                tit(i) = title(['NOISE w Envelope, f = ', num2str(freqPeak,5), ' Hz']);
                    xlim([0 1.2*(max([oddball.f standard.f noise.f]))])
                    
            i = 16;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(freqVectorWithSOA, standard.waveFFTWithSOA, 'r');
                [v,ind] = max(standard.waveFFTWithSOA);
                freqPeak = freqVectorWithSOA(ind);
                tit(i) = title(['STANDARD w SOA, f = ', num2str(freqPeak,5), ' Hz']);
                    xlim([0 1.2*(max([oddball.f standard.f noise.f]))])

            i = 17;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(freqVectorWithSOA, oddball.waveFFTWithSOA, 'b');
                [v,ind] = max(oddball.waveFFTWithSOA);
                freqPeak = freqVectorWithSOA(ind);
                tit(i) = title(['ODDBALL w SOA, f = ', num2str(freqPeak,5), ' Hz']);
                    xlim([0 1.2*(max([oddball.f standard.f noise.f]))])
            
            i = 18;
            sp(i) = subplot(rows, cols, i);
                p(i) = plot(freqVectorWithSOA, noise.waveFFTWithSOA, 'b');
                [v,ind] = max(noise.waveFFTWithSOA);
                freqPeak = freqVectorWithSOA(ind);
                tit(i) = title(['NOISE w SOA, f = ', num2str(freqPeak,5), ' Hz']);
                    xlim([0 1.2*(max([oddball.f standard.f noise.f]))])
                    
                    
            % Style plots
            set(sp, 'FontName', style.fontName, 'FontSize', style.fontBaseSize)
            set(tit, 'FontName', style.fontName, 'FontSize', style.fontBaseSize, 'FontWeight', 'bold')
            
            % Autosave the figure
            if style.imgOutautoSavePlot == 1
                fileNameOut = ['auditoryStimuli_oddballParadigm', '.png'];
                try
                    export_fig(fileNameOut, style.imgOutRes, style.imgOutAntiAlias)                
                catch err
                    err
                    warning('%s\n%s\n%s\n%s', err.identifier, 'Figure not saved most likely because you have not installed export_fig from Matlab File Exchange!', ...
                          'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"', ...
                          '"Home" -> "Environment" -> "Set Path" in Matlab2012b and newer with the updated GUI')   
                end
            end 

            
    %% SAVE WAV-files to DISK
    
        sound(oddball.wave, Fs, bits)
        pause(1.5)
        sound(oddball.waveWithEnvelope, Fs, bits)
        pause(1.5)
        sound(oddball.waveWithSOA, Fs, bits)
    
        fileName = 'standardTone.wav';
        try
            audiowrite(fileName, standard.waveWithSOA, Fs, 'BitsPerSample', bits)
        catch err
            err
            wavwrite(standard.waveWithSOA, Fs, bits, fileName)
        end
        
        fileName = 'targetTone.wav';
        try
            audiowrite(fileName, oddball.waveWithSOA, Fs, 'BitsPerSample', bits)
        catch err
            err
            wavwrite(oddball.waveWithSOA, Fs, bits, fileName)
        end
        
        fileName = 'noiseDistracter.wav';
        try
            audiowrite(fileName, noise.waveWithSOA, Fs, 'BitsPerSample', bits)
        catch err
            err
            wavwrite(noise.waveWithSOA, Fs, bits, fileName)
        end
    
            % debug on command window        
            standard
            oddball
            noise
        

            
    %% Subfunctions                    
                
    function [freqVector, ampSpectrum] = computeSingleSidedAmplitudeSpectrum(timeVector, wave, Fs)
        
        T = 1/Fs;                     % Sample time
        L = length(wave);    % Number of samples
        t = timeVector;      % Time vector        
        NFFT = 2^nextpow2(L); % Next power of 2 from length of y
        
        Y = fft(wave, NFFT) / L;
        ampSpectrum = 2*abs(Y(1:NFFT/2+1)); % single-sided amplitude spectrum
        
        freqVector = Fs/2 * linspace(0,1,NFFT/2+1);
    
    function wave = generateWhiteNoiseWithBandwidthLimitation(timeVector, fLimits, Fs)
        
        % see for example:
        % "How generating band limited white noise with matlab"
        % http://www.mathworks.com/matlabcentral/answers/80714
        
        % Generate the white noise
        % see, http://www.mathworks.com/matlabcentral/newsreader/view_thread/28239
        meanx = 0;
        stdx = 1;
        noise = meanx + stdx*randn(length(timeVector),1);        
        
            
        
        % Limit the bandwidth
        
            % design the filter
            fn = 0.5*Fs;
            fd = fLimits(1);
            fu = fLimits(2);
            fOrder = 96;
            [B1,A1] = fir1(fOrder,[fd fu]/fn);
            
            %fvtool(B1,A1,'Fs',Fs) % visualize the filter response              
            
                % Frank et al. (2012) did not specify how actually the
                % noise was bandwidth-limited (the filter order, filter
                % type), only the frequency range was given!

            % actually filter the noise
            noise_limited = filter(B1,A1,noise);
            
                % normalize to unity
                noise_limited = noise_limited / max(noise_limited);
        
        % Plot
        
            %{
            figure
        
            % plot the waveforms
            s(1) = subplot(2,2,1);
            p(1) = plot(noise);

            s(2) = subplot(2,2,2);
            p(2) = plot(noise_limited);
            
            % Plot the FFT
            
            [~, noiseFFT] = computeSingleSidedAmplitudeSpectrum(timeVector, noise, Fs);            
            [freqVector, noise_limitedFFT] = computeSingleSidedAmplitudeSpectrum(timeVector, noise_limited, Fs);
            
            s(3) = subplot(2,2,3);
            p(3) = plot(freqVector, noiseFFT);

            s(4) = subplot(2,2,4);
            p(4) = plot(freqVector, noise_limitedFFT);
            
            set(s(3:4), 'XLim', [0.2*fd 1.5*fu])
            %}
        
        % output 
        wave = noise_limited;       
       
        
        