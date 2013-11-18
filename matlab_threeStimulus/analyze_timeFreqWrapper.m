function [realCoefs, imagCoefs, timep, freq, isNaN] = analyze_timeFreqWrapper(matrixIn, parameters, epochIndex, scales, points, timeVectorIn, erpType, handles)
        
    % For the time-frequency analysis we can use the fastwavelet.m
    % provided by the ERPWaveLab, other option would be to use "cwt"
    % from Matlab's Wavelet Toolbox if you have a license for it
    % see e.g. http://www.bsp.brain.riken.jp/~phan/nfea/download.html     
    
    
    debugMatFileName = 'tempTimeFreqWrapper.mat';
    if nargin == 0
        %load('debugPath.mat')
        %load(fullfile(path.debugMATs, debugMatFileName))
        load(debugMatFileName)
        close all
    else
        save(debugMatFileName)
        %{
        if handles.flags.saveDebugMATs == 1
            path = handles.path;
            save('debugPath.mat', 'path')
            save(fullfile(path.debugMATs, debugMatFileName))            
        end
        %}
    end 
    
    % handles
    
    disp(['        .. analyzing epochs'])
    [noOfSamples, noOfChannels, noOfEpochs] = size(matrixIn);
    
    parameters.timeFreq.windowEpochs = 0;
    parameters.timeFreq.plotEpochs = 1;
    
    ITPC=0;
    ITLC=0;
    ITLCN=0;
    ERSP=0;
    avWT=0;
    WTav=0;
    avWTi=0;
    WTavi=0;
    
    % matrixIn
    %    - rows - number of samples per epoch
    %    - cols - number of channels
    
    %% Define parameters
    
        % bandwidth parameter (σ = 1) for example used in the tutorial:
        % http://www.erpwavelab.org/tutorial/index_files/Page496.htm        
        timeRes = 1 / parameters.EEG.srate;
        
        % signal would be noOfSamples (rows) x noOfChannels (cols)
        signal = matrixIn;

            % window the signal epoch
            r = parameters.timeFreq.tukeyWindowR;
            window = (tukeywin(noOfSamples,r)); % Tukey window with r=0.10 is 10% cosine window
            
            if parameters.timeFreq.windowEpochs == 1
                disp(['          ..windowing the epochs with a ', num2str(r), '% Tukey (Cosine) window'])
                for ep = 1 : noOfEpochs
                    for ch = 1 : noOfChannels
                        signal(:,ch,ep) = signal(:,ch,ep) .* window; 
                        % check if this is correct with wavelets actually?
                    end
                end                           
            end
        
            if parameters.timeFreq.plotEpochs == 1
                
                rows = 3;
                cols = noOfChannels;
                
                scrsz = get(0,'ScreenSize'); % get screen size for plotting
                fig = figure('Color', 'w');
                    set(fig, 'Position', [0.05*scrsz(3) 0.20*scrsz(4) 0.92*scrsz(3) 0.62*scrsz(4)])
                
                for ch = 1 : noOfChannels
                    
                    sp(1,ch) = subplot(rows,cols,ch);
                    hold on                 

                        trials = linspace(1, noOfEpochs, noOfEpochs);
                        epochsPerChannel = (squeeze(signal(:,ch,:)))';                    
                        
                        [T, TR] = meshgrid(timeVectorIn, trials);     
                        s = surf(T, TR, epochsPerChannel);                                  
           
                        set(s, 'EdgeColor', 'none')             
                        xlabel('Time [s]'); ylabel('Epoch #');
                        title(['Ch: ', num2str(ch)])
                        % colorbar
                        
                        pos2D(1,ch,:) = get(gca, 'Position');
                        
                        %ans =
                        %0.3361    0.7093    0.1566    0.2157

                    sp(2,ch) = subplot(rows,cols,cols+ch);
                    
                        % whos
                        epochAverage = nanmean(epochsPerChannel,1);
                        plot(timeVectorIn, epochAverage);
                        xlabel('Time [s]'); ylabel('\muV');
                        title(['Averaged waveform'])
                        drawnow

                        pos2D(2,ch,:) = get(gca, 'Position');
                        
                end
                
                
            end
            
        wname = 'cmor'; % the name of the Wavelet, continuous, Morlet

        % Scales is the parameter b, 
        fb = parameters.timeFreq.bandwidthParameter;        
        resol = 12; % 12 by default        
            
        parameters.timeFreq.timeResolutionDivider = 32;
        parameters.timeFreq.freqCutOff = [1 30];
        parameters.timeFreq.freqDownsampleFactor = 1;        

        timep = (linspace(min(timeVectorIn), max(timeVectorIn), length(points)))';

        freqIndicesFound = 0;

        %power = zeros(noOfChannels, length(scales), length(points), noOfEpochs);
        for ch = 1 : noOfChannels
            
            perChannelEpochs = squeeze(signal(:,ch,:));            
            [noOfSamples, noOfEpochs] = size(perChannelEpochs);
    
            artifactFreeEpochsFound = 0;

            for ep = 1 : noOfEpochs
                
                if ch == 1
                    isNaN(ep) = logical(sum(isnan(perChannelEpochs(:,ep))));
                    % use same indices for all the channels
                end              

                % calculate TF only for clean epochs
                if ~isNaN(ep)
                    artifactFreeEpochsFound = artifactFreeEpochsFound + 1;
                    [powerRaw, freq, scale, Cdelta, n, dj, dt, variance, coiRaw] = analyze_waveletWrapper(perChannelEpochs(:,ep), parameters.EEG.srate, parameters.timeFreq.timeResolutionDivider);
                    % isNaN = logical((squeeze(sum(squeeze((sum(isnan(squeeze(powerRaw)),2))),1))))
                    
                    % we can reduce the memory requirements slightly by cutting frequencies here
                    if freqIndicesFound == 0
                        indicesTrim = analyze_getTFtrimIndices(freq, powerRaw, coiRaw, parameters);                                             
                        freqIndicesFound = 1;
                    end

                    powerTrim = powerRaw(indicesTrim(1):indicesTrim(2),:); % trim the frequency dimension
                    freq = freq(indicesTrim(1):indicesTrim(2));

                    %disp([min(freq) max(freq)])
                    %disp([min(timep) max(timep)])
                    
                    noOfNaNs_beforeInterpolation = sum(sum(isnan(powerTrim)));                    
                    powerDownsampled = interp2(timeVectorIn, freq, powerTrim,timep,freq);                      
                    noOfNaNs_afterInterpolation = sum(sum(isnan(powerDownsampled)));
                    
                    if noOfNaNs_afterInterpolation > noOfNaNs_beforeInterpolation
                        warning('Number of NaNs increased after interpolation!')
                    elseif noOfNaNs_afterInterpolation ~= 0
                        warning('Some NaNs after interpolation!')
                    else
                        power(ch,artifactFreeEpochsFound,:,:) = powerDownsampled;
                    end

                    %coi(ch,ep,:,:)
                    %pause

                    % write out the results for better human readability of the
                    % code
                    [noOfChas, noOfEpchs, noOfScales, noOfTimePoint] = size(power);
                    powerPerChannel = squeeze(power(ch,:,:,:)); % get rid of the singleton channel dimension
                    [noOfEpchs, noOfScales, noOfTimePoint] = size(powerPerChannel);

                    averPowerPerChannel = squeeze(nanmean(powerPerChannel(:,:,:),1));
                    [noOfScales, noOfSamples] = size(averPowerPerChannel);

                else

                    % don't do anything for the artifacted epoch

                end                      
                
            end

            averageOfTheChannel(ch,:,:) = averPowerPerChannel;

        end              

        disp(['timeDiv: ', num2str(parameters.timeFreq.timeResolutionDivider), ', min f: ', num2str(min(freq)), ', max f: ', num2str(max(freq)), ', freqRes: ', num2str(freq(2)-freq(1)), ' Hz'])       

        isNaN = logical(isNaN);        
        [noOfChannels, noOfScales, noOfTimePoints, noOfEpochs] = size(power);                   
        
        %timep = (linspace(min(timeVectorIn), max(timeVectorIn), length(timeVectorIn)));        
        [T, F] = meshgrid(timep, freq);
                
        whos

        % Plot the TF
        if parameters.timeFreq.plotEpochs == 1
            
            for ch = 1 : noOfChannels

                sp(3,ch) = subplot(rows,cols,(cols*2)+ch);
                contourLevels = 64;            
                    realCoefs(ch, :, :) = real(squeeze(averageOfTheChannel(ch,:,:)));
                    imagCoefs(ch, :, :) = imag(squeeze(averageOfTheChannel(ch,:,:)));                    
                    %debug_plotTF(scales, timep, freq, realCoefs, imagCoefs, T, F, handles) 
                    
                    handle_contours = contourf(T, F, log(squeeze(realCoefs(ch, :, :))), 64, 'EdgeColor', 'none');
                    xlabel('Time [ms]'); ylabel('Frequency [Hz]');
                    title(['n = ', num2str(sum(isNaN == 0))])
                    ylim([1 30])
                    colorbar

                    pos2D(3,ch,:) = get(gca, 'Position');
                    %pos2D(1,ch,:)
                    %pos2D(2,ch,:)
                    %set(sp(1), 'Position', [pos2D(1,ch,1) pos2D(1,ch,2) pos2D(3,ch,3) pos2D(1,ch,4)])
                    %set(sp(2), 'Position', [pos2D(2,ch,1) pos2D(2,ch,2) pos2D(3,ch,3) pos2D(2,ch,4)])
                    %0.7484    0.1100    0.1146    0.2157                    

            end
            
            try
                if handles.figureOut.ON == 1    
                    drawnow
                    dateStr = plot_getDateString(); % get current date as string
                    %cd(path.outputFigures)            
                    fileNameOut = ['timeFreqDebug_',strrep(handles.inputFile, '.pdf', ''), '_', erpType, dateStr];

                    disp([' ... saving figure to disk (', fileNameOut, '.png]'])
                    export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig(1))
                    %cd(path.code)
                end
            catch err
                err
            end
            
        end

        %% CHECK OTHER DERIVED measures from tfanalysis.m of ERPWAVELAB


    function indicesTrim = analyze_getTFtrimIndices(f, powerRaw, coiRaw, parameters)

        % low cut
        indexTemp = find(f == parameters.timeFreq.freqCutOff(1));           
        if isempty(~indexTemp)      
            freqDiff = (f - parameters.timeFreq.freqCutOff(1));
            freqDiff(freqDiff > 0) = NaN;
            [val1, fixed_ind1] = min(freqDiff);
            disp(['    The closest match is = ', num2str(f(fixed_ind1)), ' Hz (index = ', num2str(fixed_ind1), ')'])
            indexTemp = fixed_ind1;
        end
        indicesTrim(1) = indexTemp; % only one value should be found

        % hi cut
        indexTemp = find(f == parameters.timeFreq.freqCutOff(2));           
        if isempty(~indexTemp)      
            freqDiff = (f - parameters.timeFreq.freqCutOff(2));
            % get the minimum of values above zero to ensure that the
            % frequency picked is higher than the high cutoff
            freqDiff(freqDiff < 0) = NaN;
            [val2, fixed_ind2] = min(freqDiff);
            disp(['    The closest match is = ', num2str(f(fixed_ind2)), ' Hz (index = ', num2str(fixed_ind2), ')'])
            indexTemp = fixed_ind2;
        end
        indicesTrim(2) = indexTemp; % only one value should be found

