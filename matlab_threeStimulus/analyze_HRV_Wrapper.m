function [heart, hp4] = analyze_HRV_Wrapper(heart, hpRaw, thpRaw, hp, thp, rrTimes, rrPeakInterval, rrPeakAmplitude, heartrateSampleRate, parameters, handles)

    debugMatFileName = 'tempHRV.mat';
    if nargin == 0
        load('debugPath.mat')
        load(fullfile(path.debugMATs, debugMatFileName))
        close all
    else
        if handles.flags.saveDebugMATs == 1
            path = handles.path;
            save('debugPath.mat', 'path')
            save(fullfile(path.debugMATs, debugMatFileName))            
        end
    end         

    debugON = 1;   

    % get stats [units in milliseconds, x1000]
    avIBI = mean(hp*1000);
    maxIBI = max(hp*1000);
    minIBI = min(hp*1000);
    RMS = RMSSD(hp*1000);
    SDNN = std(hp*1000);

    hp = hp*1000;       

    if ~strcmp(parameters.heart.detrendMethod, 'ougp')
    %% CLASSICAL INTERPOLATION METHOD
       % spline interpolation
        disp(['              - Spline interpolation to ', num2str(parameters.heart.rrTimeSeriesUpsampleRate), ' Hz'])
        auxtime = thp(1) : 1/parameters.heart.rrTimeSeriesUpsampleRate : thp(end);
        hp2 = interp1(thp, hp, auxtime, 'spline')';    
        % hp2 = pre_nyquistInterpolation(hp, length(auxtime), 1,0); % Nyquist interpolation  

            % See Fisher et al. (2012)  for the discussion of interpolation,
            % and an improved method to avoid this step

                % Fisher AC, Eleuteri A, Groves D, Dewhurst CJ. 2012. 
                % The Ornstein–Uhlenbeck third-order Gaussian process (OUGP) applied directly to the un-resampled heart rate variability (HRV) tachogram for detrending and low-pass filtering. 
                % Med Biol Eng Comput 50:737–742. 
                % http://dx.doi.org/10.1007/s11517-012-0928-2
                % http://clinengnhs.liv.ac.uk/OUGP.htm        
    end

    % detrend hp
    if strcmp(parameters.heart.detrendMethod, 'constant')
        hp3 = detrend(hp2, parameters.heart.detrendMethod);            
        useOUGP = 0;
    elseif strcmp(parameters.heart.detrendMethod, 'linear')
        hp3 = detrend(hp2, parameters.heart.detrendMethod);
        useOUGP = 0;
    elseif strcmp(parameters.heart.detrendMethod, 'ougp')
        % Ornstein-Uhlenbeck 3rd order Gaussian Process method
        % http://dx.doi.org/10.1007/s11517-012-0928-2

        % with OUGP we are going to bandpass filter the input for the
        % Frequencies of Interest of HRV output (thus the workflow slightly
        % differs from the Kardia implementation and it is easier to
        % maintain the old switches here and use the boolean flag below to
        % calculate with OUGP method
        useOUGP = 1;
        disp(['              - Ornstein–Uhlenbeck third-order Gaussian process (OUGP) filtering'])
    else
        detrendMethod = parameters.heart.detrendMethod;
        error('what detrending for heart spectral analysis you wanted? Typo?')
    end

    if useOUGP == 0

        % Filter method
        if strcmp(parameters.heart.filterWindow, 'hanning')
            wdw=hanning(length(hp3));        
        elseif strcmp(parameters.heart.filterWindow, 'hamming')
            wdw=hamming(length(hp3));
        elseif strcmp(parameters.heart.filterWindow, 'blackman')        
            wdw=blackman(length(hp3));
        elseif strcmp(parameters.heart.filterWindow, 'bartlett')    
            wdw=bartlett(length(hp3));
        else
            wdw = parameters.heart.filterWindow
            error('what window for heart spectral analysis you wanted? Typo?')
        end


        hp4=hp3.*wdw;

        % Calculate FFT points
        % the 'very low frequency' band is given at 0.04 Hz, thus the NFFT
        % should be at least 25s (1 / 0.04) to get that kind of spectral
        % resolution


        if strcmp(parameters.heart.noOfFFTPoints, '512')
            N=512;        
        elseif strcmp(parameters.heart.noOfFFTPoints, '1024')
            N=1024;
        elseif strcmp(parameters.heart.noOfFFTPoints, 'auto')
            N = 2^nextpow2(length(hp4));
        else
            N = 2^nextpow2(length(hp4))
            error('how many points for FFT of heart spectral analysis you wanted? Typo?')
        end      

        % Original KARDIA Implementation of frequency analysis
        fs = heartrateSampleRate;
        if strcmp(parameters.heart.spectralAlgorithm, 'FFT')
            cw = (1/N) * sum(wdw.^2); % Coefficient to remove window effect
            PSD = (abs(fft(hp4,N)).^2)/(N*fs*cw);        
            F = (0:fs/N:fs-fs/N)';        
            PSD = 2*PSD(1:ceil(length(PSD)/2));        
            F = F(1:ceil(length(F)/2));

            % also for non-interpolated
            if debugON == 1
                % have to implement later
            end

        elseif strcmp(parameters.heart.spectralAlgorithm, 'AR model')
            [A, variance] = arburg(hp4,ARorder);
            [H,F] = freqz(sqrt(variance),A,N/2,fs);
            cw = (1/length(hp4)) * sum(wdw.^2); % Coefficient to remove the window effect
            PSD = 2*(abs(H).^2)/(fs*cw);         

            % also for non-interpolated
            if debugON == 1
                % have to implement later
            end

        else
            spectralAlgorithm = parameters.heart.spectralAlgorithm
            error('What spectral analysis you wanted for heart rate? typo?')
        end

        % get power in different bands    
        ULF = spPCRower(F, PSD, 'ulf', 'kardia'); % added by Petteri to match the bands from OUGP paper
        ULFStar = spPCRower(F, PSD, 'ulfStar', 'kardia'); % added by Petteri to match the bands from OUGP paper
        VLF = spPCRower(F, PSD, 'vlf', 'kardia');
        LF = spPCRower(F, PSD, 'lf', 'kardia');
        HF = spPCRower(F, PSD, 'hf', 'kardia');
        nlf = spPCRower(F, PSD, 'nlf', 'kardia');
        nhf = spPCRower(F, PSD, 'nhf', 'kardia');

    
    elseif useOUGP == 1                       
    %% OUGP based non-interpolation method
                
        % BIN-by-BIN
        if parameters.heart.ougpBinByBin == 1
            freqBinNames = fieldnames(parameters.heart.freqBins);
            noOfFreqBins = length(freqBinNames);

            for bin =  1 : noOfFreqBins

                bandBassFreqs = parameters.heart.freqBins.(freqBinNames{bin});
                %thp'
                %hp'

                % Detrend and band-pass filter
                [smooth.(freqBinNames{bin}), nf.(freqBinNames{bin})] = ougp_gpousmooth2(thp', hp', bandBassFreqs, 'bp');

                % Calculate the PSD now using Lomb-Scargle instead of FFT or AR
                % as done with KARDIA
                F.(freqBinNames{bin}) = bandBassFreqs(1) : parameters.heart.freqResolution : bandBassFreqs(2);
                [Pr, fr, A, z0, A0, ofac] = ougp_flspw(thp, smooth.(freqBinNames{bin}), F.(freqBinNames{bin}));
                PSD.(freqBinNames{bin}) = A;
                [bin size(F.(freqBinNames{bin})) size(PSD.(freqBinNames{bin})) size(fr{1})]

                % have to combine the different bins somehow still so that
                % the band power estimation won't crash
                
            end
        else            
            try
                [hp4, nyqFreq] = ougp_gpousmooth2(thp, hp, parameters.heart.freqRange, 'bp');
                    % if you get only NaNs out, check the freqRange used
            catch err
                % err
                % whos
                disp('                .. transpose thp and hp to avoid crash')
                [hp4, nyqFreq] = ougp_gpousmooth2(thp', hp', parameters.heart.freqRange, 'bp');
            end
            F = parameters.heart.freqRange(1) : parameters.heart.freqResolution : parameters.heart.freqRange(2);
            
            % for the input
            if debugON == 1
                [Pr, fr, PSD_In, z0, A0, ofac] = ougp_flspw(thp, hp, F);
            end
            
            % for the filtered
            [Pr, fr, A, z0, A0, ofac] = ougp_flspw(thp, hp4, F);
            PSD = A;          
            
        end
        
        F = F(1:end-1); % check if correct, that the last value is clipped indeed, and not the first
        
        % get power in different bands    
        ULF = spPCRower(F, PSD, 'ULF', parameters, 'OUGP'); % added by Petteri to match the bands from OUGP paper
        ULFStar = spPCRower(F, PSD, 'ULFStar', parameters, 'OUGP'); % added by Petteri to match the bands from OUGP paper
        VLF = spPCRower(F, PSD, 'VLF', parameters, 'OUGP');
        LF = spPCRower(F, PSD, 'LF', parameters, 'OUGP');
        HF = spPCRower(F, PSD, 'HF', parameters, 'OUGP');
        TOTAL = spPCRower(F, PSD, 'TOTAL', parameters, 'OUGP');
        
    else
        
    end

    
    
    
    %% DEBUG
    if debugON == 1
        scrsz = handles.style.scrsz;
        fig = figure('Color', 'w', 'Name', 'HRV');
            set(fig, 'Position', [0.25*scrsz(3) 0.05*scrsz(4) 0.45*scrsz(3) 0.90*scrsz(4)])
            rows = 3; cols = 1;
        
        % Original time series
        ind = 1; sp(ind) = subplot(rows, cols, ind);
            hold on
            plot(thpRaw, 1000*hpRaw, 'r')
            plot(thp, hp, 'b'); 
            hold off
                legend('Outlier Rejected', 'Original Raw', 'Location', 'Best')
                legend('boxoff')
            title('Original R-R series')
            xlim([0 max(thp)])
            lab(ind,1) = xlabel('Time [s]');
            lab(ind,2) = ylabel('RR Interval [ms]');
        
        % Filtered time series    
        ind = 2; sp(ind) = subplot(rows, cols, ind);
        
            if useOUGP == 1       
                
                if parameters.heart.ougpBinByBin == 1
                    filteredRR = zeros(noOfFreqBins, length(thp));
                    for bin =  1 : noOfFreqBins
                        % make matrix
                        filteredRR(bin,:) = smooth.(freqBinNames{bin});
                    end
                    plot(thp', filteredRR); title('Detrended and bandpass filtered R-R series')
                    leg(1) = legend(freqBinNames);
                    set(leg, 'Position',[0.856150793650792 0.558967391304347 0.131613756613757 0.116847826086957])
                    legend('boxoff')
                else
                    plot(thp', hp4); title(['Detrended and bandpass filtered R-R series, f_N_y_q = ', num2str(nyqFreq, '%1.3f'), ' Hz'])
                end
                
            else
                plot(auxtime, hp2); title(['R-R series upsampled to ', num2str(parameters.heart.rrTimeSeriesUpsampleRate), ' Hz'])
            end
            xlim([0 max(thp)])
            lab(ind,1) = xlabel('Time [s]');
            lab(ind,2) = ylabel('RR Interval [\Deltams]');
            
        % Power Spectral Density (PSD)
        ind = 3; sp(ind) = subplot(rows, cols, ind);
        
            if useOUGP == 1
                
                % Bin-by-bin
                if parameters.heart.ougpBinByBin == 1
                    hold on
                    for bin =  1 : noOfFreqBins
                        % make matrix                    
                        % OUGP_PSD(bin,:) = PSD.(freqBinNames{bin});
                        x = F.(freqBinNames{bin})(1:end-1);
                        y = PSD.(freqBinNames{bin});
                        plot(x, y);
                    end
                    hold off
                    leg2 = legend(freqBinNames);
                    set(leg2, 'Position',[0.848214285714283 0.235054347826087 0.131613756613757 0.116847826086957]);
                    legend('boxoff')
                else
                    %size(F)
                    %size(PSD)
                    hold on                    
                    area(F, PSD_In, 'FaceColor', 'r', 'EdgeColor', 'none');
                    plot(F, PSD, 'b', 'LineWidth', 1)
                    leg2 = legend('Filtered', 'Input');
                    %set(leg2, 'Position',[0.848214285714283 0.235054347826087 0.131613756613757 0.116847826086957]);
                    legend('boxoff')
                    
                end
                title('Lomb-Scargle PSD')
                
            else
                plot(F, PSD); xlim([0 0.5]); title('PSD')
            end
            lab(ind,1) = xlabel('Frequency [Hz]');
            lab(ind,2) = ylabel('Power [s^2Hz^{-1}]');
            
            % Auto-SAVE
            try
                if handles.figureOut.ON == 1      
                    drawnow
                    dateStr = plot_getDateString(); % get current date as string          
                    fileNameOut = sprintf('%s%s', 'debug_HRV_', strrep(handles.inputFile, '.bdf', ''), '.png');
                    export_fig(fullfile(handles.path.debugHeartOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
                    %cd(path.code)
                end
            catch err
                err
                str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                              'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
                error(str)
            end
    end
    

    % save results structure
    heart.vector.RR_t = thp;
    heart.vector.RR_timeSeriesAfterOutlier = hp;    
    heart.vector.RR_timeSeriesFilt = hp4;
    heart.vector.PSD_In = PSD_In;
    heart.vector.PSD = PSD;
    heart.vector.F = F;

    heart.scalar.TOTAL = TOTAL;
    heart.scalar.HF = HF;
    heart.scalar.LF = LF;
    heart.scalar.VLF = VLF;
    heart.scalar.LF_HF_ratio = LF / HF;
    heart.scalar.NHF = HF / (HF+LF);
    heart.scalar.NLF = LF / (HF+LF);
    heart.scalar.avIBI = avIBI;
    heart.scalar.maxIBI = maxIBI;
    heart.scalar.minIBI = minIBI;
    heart.scalar.RMSSD = RMS;
    heart.scalar.SDNN = SDNN;


