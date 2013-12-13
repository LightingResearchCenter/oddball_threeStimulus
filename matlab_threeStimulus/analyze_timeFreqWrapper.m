function [realCoefs, imagCoefs, realCoefs_SD, imagCoefs_SD, timep, freq, isNaN, allEpochs, derivedMeasures] = analyze_timeFreqWrapper(matrixIn, parameters, epochIndex, scales, points, timeVectorIn, erpType, IAF_peak, handles)
        
    % For the time-frequency analysis we can use the fastwavelet.m
    % provided by the ERPWaveLab, other option would be to use "cwt"
    % from Matlab's Wavelet Toolbox if you have a license for it
    % see e.g. http://www.bsp.brain.riken.jp/~phan/nfea/download.html         
    
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
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
    end
    
    % handles
    
    disp(['        .. analyzing epochs'])
    [noOfSamples, noOfChannels, noOfEpochs] = size(matrixIn);
    
    parameters.timeFreq.windowEpochs = 1;
    parameters.timeFreq.plotEpochs = 1;
    
    handles.style.lineGrey = [0.4 0.4 0.4];
    
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
            r = 2*r;
            window = (tukeywin(noOfSamples,r)); % Tukey window with r=0.10 is 10% cosine window
            
            if parameters.timeFreq.windowEpochs == 1
                disp(['          ..windowing the epochs with a ', num2str(100*r), '% Tukey (Cosine) window'])
                for ep = 1 : noOfEpochs
                    for ch = 1 : noOfChannels
                        signal(:,ch,ep) = signal(:,ch,ep) .* window; 
                        % check if this is correct with wavelets actually?
                    end
                end  
            else
                disp(['          ..not windowing the epochs with a ', num2str(100*r), '% Tukey (Cosine) window (for example)'])
            end
        
            %% PLOT INPUT EPOCHS AND AVERAGE WAVEFORM
            if parameters.timeFreq.plotEpochs == 1
                
                rows = 3;
                cols = noOfChannels;
                
                scrsz = get(0,'ScreenSize'); % get screen size for plotting
                fig = figure('Color', 'w', 'Name', ['TF: ', erpType]);
                    set(fig, 'Position', [0.05*scrsz(3) 0.20*scrsz(4) 0.92*scrsz(3) 0.62*scrsz(4)])
                
                for ch = 1 : noOfChannels
                
                    % EPOCHS
                    sp(1,ch) = subplot(rows,cols,ch);
                    hold on                 

                        trials = linspace(1, noOfEpochs, noOfEpochs);
                        epochsPerChannel = (squeeze(signal(:,ch,:)))';                    
                        epochLimits(1,ch,:) = [min(min(epochsPerChannel)) max(max(epochsPerChannel))];
                        
                        [T, TR] = meshgrid(timeVectorIn, trials);     
                        s(ch) = surf(T, TR, epochsPerChannel);                                        
                        caxis manual                        
                        
                        lab(ch,1,1) = xlabel('Time [s]'); 
                        lab(ch,2,1) = ylabel('Epoch #');
                        tit(ch,1) = title(['Ch: ', num2str(ch), ' (', parameters.BioSemi.chName{ch+parameters.BioSemi.chOffset}, ')']);                                                
                        
                        if ch == noOfChannels
                            cb = colorbar('peer',gca,...
                                [0.92415912031048 0.712749615975423 0.0161707632600259 0.204301075268817]);
                            set(get(cb,'title'),'String','\muV');

                        end
                        ylim([1 noOfEpochs])
                        pos2D(1,ch,:) = get(gca, 'Position');
                        
                        %ans =
                        %0.3361    0.7093    0.1566    0.2157

                    % AVERAGE WAVEFORM
                    sp(2,ch) = subplot(rows,cols,cols+ch);
                    
                        % whos
                        dim = 1;
                        epochAverage = nanmean(epochsPerChannel,dim);
                        epochSD = nanstd(epochsPerChannel,0,dim);                        
                        epochLimits(2,ch,:) = [min(epochAverage-epochSD) max(epochAverage+epochSD)];
                        
                        alpha = .4;                        
                        pSh(ch,:) = plot_errorShade(timeVectorIn, epochAverage', epochSD', alpha, [0 0.60 1]);
                        
                        lab(ch,1,2) = xlabel('Time [s]'); 
                        lab(ch,2,2) = ylabel('\muV');                        
                        tit(ch,2) = title(['Averaged waveform']);                        

                        pos2D(2,ch,:) = get(gca, 'Position');
                        
                end
                
                % STYLE
                for ch = 1 : noOfChannels
                    % fix the colorbar limits to be the same
                    try
                        caxis(sp(1,ch), [min(min(epochLimits(1,:,:))) max(max(epochLimits(1,:,:)))])
                    catch err
                        err
                        error('Probably no epochs plotted so you cannot fix the colorbar limits, handle earlier that the script "knows" that there are no artifact-free epochs')
                    end
                end
                set(s, 'EdgeColor', 'none')         
                yLims = [min(min(epochLimits(2,:,:))) max(max(epochLimits(2,:,:)))];
                set(sp(2,:), 'YLim', yLims)
                
                for ch = 1 : noOfChannels
                    axes(sp(1,ch))
                    zeroLineH(1,ch) = line([0 0], [1 noOfEpochs], 'Color', handles.style.lineGrey);                    
                    axes(sp(2,ch))
                    zeroLineH(2,ch) = line([0 0], yLims, 'Color', handles.style.lineGrey);
                    zeroLineH(3,ch) = line([min(timeVectorIn) max(timeVectorIn)], [0 0], 'Color', handles.style.lineGrey);
                    % push the average back to front
                    uistack(pSh(ch,3), 'top')
                end
                drawnow
                
            end
            
        %% WAVELET TRANSFORM
            
            wname = 'cmor'; % the name of the Wavelet, continuous, Morlet

            % Scales is the parameter b, 
            fb = parameters.timeFreq.bandwidthParameter;        
            resol = 12; % 12 by default        
            
            % parameters.timeFreq.timeResolutionDivider = 32;
            parameters.timeFreq.freqCutOff = [1 30];
            parameters.timeFreq.freqDownsampleFactor = 1;        

            timep = (linspace(min(timeVectorIn), max(timeVectorIn), length(points)))';

            freqIndicesFound = 0;

            % find artifacted epochs
            for ch = 1 : noOfChannels 
                isNaN(:,ch) = logical(sum(isnan(squeeze(signal(:,ch,:)))))';
            end

            fprintf('         - ')
            %power = zeros(noOfChannels, length(scales), length(points), noOfEpochs);
            for ch = 1 : noOfChannels

                fprintf(['ch', num2str(ch), ' '])
                perChannelEpochs = squeeze(signal(:,ch,~isNaN(:,ch)));
                [noOfSamples, noOfEpochs] = size(perChannelEpochs);

                %power
                if noOfEpochs ~= 0
                    
                    for ep = 1 : noOfEpochs
                        
                        % The wavelet transform is now computed on a
                        % single-trial basis for each epoch and the
                        % individual TF epochs are then averaged (or
                        % manipulated as wished)
                        
                        % Visual representation of why is this so, one can
                        % check out the Figure 4 of Herrmann et al. (2013)
                        % "Time–Frequency Analysis of Event-Related Potentials: A Brief Tutorial"
                        % http://dx.doi.org/10.1007/s10548-013-0327-5
                        
                        % In other words, we would not necessarily make a
                        % huge error with the time-frequency spectrum of
                        % the evoked potential, but we might lose a lot of
                        % information on the induced potential

                        %% WAVELET WRAPPER
                        [powerRaw, WT, freq, scale, Cdelta, n, dj, dt, variance, coiRaw] = analyze_waveletWrapper(perChannelEpochs(:,ep), parameters.EEG.srate, parameters.timeFreq.timeResolutionDivider);                                         
                        
                        % we can reduce the memory requirements slightly by cutting frequencies here
                        if freqIndicesFound == 0
                            indicesTrim = analyze_getTFtrimIndices(freq, [parameters.timeFreq.freqCutOff(1) parameters.timeFreq.freqCutOff(2)], [], [], parameters);                                             
                            freqIndicesFound = 1;
                        end

                        WTtrim = WT(indicesTrim(1):indicesTrim(2),:); % trim the frequency dimension
                        power = powerRaw(indicesTrim(1):indicesTrim(2),:); % trim the frequency dimension
                        freq2 = freq(indicesTrim(1):indicesTrim(2));
                        %disp([min(freq) max(freq)])
                        %disp([min(timep) max(timep)])                       
                        
                        powerDownsWithCOI = interp2(timeVectorIn, freq2, power, timep, freq2);
                        wtDownsWithCOI = interp2(timeVectorIn, freq2, WTtrim, timep, freq2);
                        
                        % Reject low frequencies under the COI:                                                
                        % a vector of N points that contains the maximum period of useful information
                        % at that particular time. Periods greater than this are subject to edge effects.                        
                        % plot_debugCOI_onWaveletTransform(timeVectorIn, freq, 1./freq, powerRaw, coiRaw, handles)
                        [power, nanMask] = analyze_rejectWT_underCOI(timeVectorIn, freq2, power, coiRaw, parameters, handles);
                        [WTtrim, nanMask] = analyze_rejectWT_underCOI(timeVectorIn, freq2, WTtrim, coiRaw, parameters, handles);
                        
                        % downsample the time resolution, keeping the frequency
                        % vector unchanged
                        % noOfNaNs_beforeInterpolation = sum(sum(isnan(WTtrim)));                    
                        
                        % now the WT and power contain NaNs as we rejected
                        % some values based on the COI                                                
                        power = interp2(timeVectorIn, freq2, power, timep, freq2);
                        WT_Downsampled = interp2(timeVectorIn, freq2, WTtrim, timep, freq2);                        
                        
                        
                        % noOfNaNs_afterInterpolation = sum(sum(isnan(WT_Downsampled)));                                                
                        % coi_downsampled = interp1(timeVectorIn, coiRaw', timep);% Cone-of-Influence
                        

                        %% NORMALIZE if needed
                        normalizeTheWaveletSpectrum = 1; % could put this to init_DefaultParameters.m
                                                         % but you basically
                                                         % always want to
                                                         % normalize?                    

                            if normalizeTheWaveletSpectrum == 1      
                                
                                debugOn = 0;
                                
                                [power, nonNormFreqIndex] = analyze_normalizeWaveletSpectrum(power, powerDownsWithCOI, timep, freq, ...
                                                    parameters.oddballTask.ERP_baselineCorrection, parameters.timeFreq.timeResolutionDivider, debugOn, parameters, handles);
                                                
                                [WT_Downsampled, nonNormFreqIndex] = analyze_normalizeWaveletSpectrum(WT_Downsampled, wtDownsWithCOI, timep, freq, ...
                                                    parameters.oddballTask.ERP_baselineCorrection, parameters.timeFreq.timeResolutionDivider, debugOn, parameters, handles);                                                
                                
                            end

                        % Check if the interpolation failed           
                        WTout{ch}(ep,:,:) = WT_Downsampled;
                        powerOut{ch}(ep,:,:) = power;                        

                        % EXTRA OUTPUTS (PARAMETERS CALCULATED)
                        %whos('WT*', 'IT*', 'ERSP', 'avWT', 'WTav', 'avWTi', 'WTavi', 'parameters', 'handles')
                        [ITPC, ITLC, ITLCN, ERSP, avWT, WTav, avWTi, WTavi] = ...
                            analyze_wavelet_derivedParametersAccum(WT_Downsampled, ITPC, ITLC, ITLCN, ERSP, avWT, WTav, avWTi, WTavi, parameters, handles);
                        
                        % write out the results for better human readability of the code                   
                        % [noOfEpochs2, noOfScales, noOfTimePoints] = size(power);
                        averPowerPerChannel = squeeze(nanmean(powerOut{ch}(:,:,:),1));
                        stdOfPowerPerChannel = squeeze(nanstd(powerOut{ch}(:,:,:),1));
                        % [noOfScales, noOfTimePoints] = size(averPowerPerChannel);
                        % disp([noOfScales noOfTimePoint])


                    end

                    [derivedMeasures.ITPC{ch}, derivedMeasures.ITLC{ch}, derivedMeasures.ITLCN{ch}, derivedMeasures.ERSP{ch}, derivedMeasures.avWT{ch}, derivedMeasures.WTav{ch}, derivedMeasures.Induced{ch}] = ...
                        analyze_wavelet_derivedParameters(ITPC, ITLC, ITLCN, ERSP, avWT, WTav, avWTi, WTavi, noOfEpochs, parameters, handles);
                    
                    % now the different channels might have different amount of
                    % epochs so we cannot put all the values to a single matrix, so
                    % we use a cell
                    %averPowerPerChannel = nanmean(WTout(:,:,:),1);
                    %stdOfPowerPerChannel = nanstd(WTout(:,:,:),1);
                    averageOfTheChannel{ch} = averPowerPerChannel;
                    stdOfTheChannel{ch} = stdOfPowerPerChannel;

                else                    
                    
                    powerOut{ch} = NaN;
                    WTout{ch} = NaN;
                    
                    % if no valid epochs are found, something has to be done to
                    % avoid the script from crashing
                    rowsNaN = 1; colsNaN = 1;
                    [WTout{ch}, powerOut{ch}, averageOfTheChannel{ch}, stdOfTheChannel{ch}, derivedMeasures] = analyze_fillWaveletOutputWithNaNs(rowsNaN,colsNaN,ch,handles);                    
                    % disp(['                 .. no valid epochs found, NaN-filled'])
                    
                end
            end
            fprintf('\n')
            disp(['            - timeDiv: ', num2str(parameters.timeFreq.timeResolutionDivider), ' (timeRes = ', num2str(1000*(timep(2)-timep(1))), ' ms), min f: ', num2str(min(freq)), ', max f: ', num2str(max(freq)), ', freqRes: ', num2str(freq(2)-freq(1)), ' Hz'])       

            
    %% ASSIGN TO BE RETURNED    
    allEpochs.WT = WTout;
    allEpochs.power = powerOut;            
            
            
            
    %% PLOT THE TIME-FREQUENCY 
    
        isNaN = logical(isNaN);        
        [noOfScales, noOfTimePoints, noOfEpochs] = size(WTout{1});                   
        noOfChannels = length(averageOfTheChannel);
        [T, F] = meshgrid(timep, freq2);  
                   
        for ch = 1 : noOfChannels

            %% Output to the calling function
            realCoefs{ch}(:,:) = real(squeeze(averageOfTheChannel{ch}(:,:)));
                % size(realCoefs{ch})
            realCoefs_SD{ch}(:,:) = real(squeeze(stdOfTheChannel{ch}(:,:)));
            imagCoefs{ch}(:,:) = imag(squeeze(averageOfTheChannel{ch}(:,:)));                    
            imagCoefs_SD{ch}(:,:) = imag(squeeze(stdOfTheChannel{ch}(:,:)));
                
            %% Plot the TF
            if parameters.timeFreq.plotEpochs == 1
                
                sp(3,ch) = subplot(rows,cols,(cols*2)+ch);
                contourLevels = 256;
                Z = averageOfTheChannel{ch};
                
                %debug_plotTF(scales, timep, freq, realCoefs, imagCoefs, T, F, handles) 
                try
                    [c, handle_contours{ch}] = contourf(T, F, Z, contourLevels, 'EdgeColor', 'none'); 
                catch err
                    % goes here when no valid epochs are found
                    % err
                end
                
                epochLimits(3,ch,:) = [min(min(Z)) max(max(Z))];
                
                lab(ch,1,3) = xlabel('Time [ms]'); 
                lab(ch,2,3) = ylabel('Frequency [Hz]');                    
                tit(ch,3) = title(['n = ', num2str(sum(isNaN(:,ch) == 0)), ', IAF = ', num2str(IAF_peak), ' Hz']);
                ylim([1 30])

                % annotate the individual alpha range                    
                alphaLimits = IAF_peak + parameters.powerAnalysis.eegBins.freqs{1}; % 1 is now hard-coded, see init_defaultParameters.m to make this more robust
                
                    % Now we are annotating the alpha based on the IAF of
                    % each individual on each session thus "normalizing"
                    % frequency variations across subjects, additionally
                    % you could try to time-lock the epochs to the button
                    % press as there is possibly some motor-contamination
                    % of the epochs, and some of the jitter could be
                    % explained by variations in reaction time.
                    
                    % See for example:
                    % Makeig S, Onton J. 2009. ERP Features and EEG Dynamics: An ICA Perspective.
                    % http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.160.4983


                hold on
                lineHandle(ch,1) = line([min(timep) max(timep)], [alphaLimits(1) alphaLimits(1)]);
                lineHandle(ch,2) = line([min(timep) max(timep)], [alphaLimits(2) alphaLimits(2)]);
                hold off

                if ch == noOfChannels
                    cb = colorbar('peer',gca,...
                        [0.929333764553688 0.104454685099846 0.0161707632600259 0.205837173579109]);
                    set(get(cb,'title'),'String','%');

                end

                pos2D(3,ch,:) = get(gca, 'Position');
                %pos2D(1,ch,:)
                %pos2D(2,ch,:)
                %set(sp(1), 'Position', [pos2D(1,ch,1) pos2D(1,ch,2) pos2D(3,ch,3) pos2D(1,ch,4)])
                %set(sp(2), 'Position', [pos2D(2,ch,1) pos2D(2,ch,2) pos2D(3,ch,3) pos2D(2,ch,4)])
                %0.7484    0.1100    0.1146    0.2157                    

            end                     
            
        end
        
        % STYLE
        set(lineHandle, 'Color', 'w', 'LineStyle', '-')
        set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2) 
        set(sp, 'XLim', [min(timep) max(timep)]) 
        set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold') 
        set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2, 'FontWeight', 'bold')            
        
            for ch = 1 : noOfChannels
                %caxis(sp(3,ch), [min(min(epochLimits(3,:,:))) max(max(epochLimits(3,:,:)))]);
                caxis(sp(3,ch), [0 600]);
                %set(sp(3,ch), 'ZScale', 'log')
            end
            
        % AUTO-SAVE FIGURE
        try
            if handles.figureOut.debugON == 1 
                drawnow
                dateStr = plot_getDateString(); % get current date as string
                %cd(path.outputFigures)            
                fileNameOut = ['timeFreqDebug_',strrep(handles.inputFile, '.bdf', ''), '_', erpType, dateStr];

                disp(['         ... saving figure to disk (', fileNameOut, '.png]'])
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig(1))
                %cd(path.code)
            end
        catch err
            err
        end           
        


% Fills the output with NaN values in case where no epochs are found for
% the channel (e.g. when the electrode came off during session and all data
% is garbage)
function [WTout, powerOut, averageOfTheChannel, stdOfTheChannel, derivedMeasures] = analyze_fillWaveletOutputWithNaNs(rows,cols,ch,handles)

    % one could fill out the values to have the same dimensions as the
    % channels with some epohcs, but maybe it is less memory-demanding just
    % to use NaN value, and handle the possible dimension mismatch in the
    % batch part of the analysis (but still you have the rows, cols as
    % input arguments if you want to do that)
    
    WTout = zeros(rows,cols) * NaN;
    powerOut = zeros(rows,cols) * NaN;
    averageOfTheChannel = zeros(rows,cols) * NaN;
    stdOfTheChannel = zeros(rows,cols) * NaN;

    derivedMeasures.ITPC{ch} = zeros(rows,cols) * NaN;
    derivedMeasures.ITLC{ch} = zeros(rows,cols) * NaN;
    derivedMeasures.ITLCN{ch} = zeros(rows,cols) * NaN;
    derivedMeasures.ERSP{ch} = zeros(rows,cols) * NaN;
    derivedMeasures.avWT{ch} = zeros(rows,cols) * NaN;
    derivedMeasures.WTav{ch} = zeros(rows,cols) * NaN;
    derivedMeasures.Induced{ch} = zeros(rows,cols) * NaN;
    