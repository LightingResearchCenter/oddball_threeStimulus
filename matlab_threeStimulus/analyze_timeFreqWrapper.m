function [realCoefs, imagCoefs, realCoefs_SD, imagCoefs_SD, timep, freq, isNaN] = analyze_timeFreqWrapper(matrixIn, parameters, epochIndex, scales, points, timeVectorIn, erpType, IAF_peak, handles)
        
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
    
    parameters.timeFreq.windowEpochs = 0;
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
        
            %% PLOT INPUT EPOCHS AND AVERAGE WAVEFORM
            if parameters.timeFreq.plotEpochs == 1
                
                rows = 3;
                cols = noOfChannels;
                
                scrsz = get(0,'ScreenSize'); % get screen size for plotting
                fig = figure('Color', 'w');
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
            
            parameters.timeFreq.timeResolutionDivider = 32;
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
                for ep = 1 : noOfEpochs

                    %% WAVELET WRAPPER
                    [power, WT, freq, scale, Cdelta, n, dj, dt, variance, coiRaw] = analyze_waveletWrapper(perChannelEpochs(:,ep), parameters.EEG.srate, parameters.timeFreq.timeResolutionDivider);                                         

                    % we can reduce the memory requirements slightly by cutting frequencies here
                    if freqIndicesFound == 0
                        indicesTrim = analyze_getTFtrimIndices(freq, [parameters.timeFreq.freqCutOff(1) parameters.timeFreq.freqCutOff(2)], [], [], parameters);                                             
                        freqIndicesFound = 1;
                    end

                    WTtrim = WT(indicesTrim(1):indicesTrim(2),:); % trim the frequency dimension
                    power = power(indicesTrim(1):indicesTrim(2),:); % trim the frequency dimension
                    freq2 = freq(indicesTrim(1):indicesTrim(2));
                    %disp([min(freq) max(freq)])
                    %disp([min(timep) max(timep)])

                    % downsample the time resolution, keeping the frequency
                    % vector unchanged
                    noOfNaNs_beforeInterpolation = sum(sum(isnan(WTtrim)));                    
                    WT_Downsampled = interp2(timeVectorIn, freq2, WTtrim, timep, freq2);                      
                    power = interp2(timeVectorIn, freq2, power, timep, freq2);                      
                    noOfNaNs_afterInterpolation = sum(sum(isnan(WT_Downsampled)));

                    %% NORMALIZE if needed
                    normalizeTheWaveletSpectrum = 1; % could put this to init_DefaultParameters.m
                                                     % but you basically
                                                     % always want to
                                                     % normalize?                    

                        if normalizeTheWaveletSpectrum == 1 
                            WT_Downsampled = analyze_normalizeWaveletSpectrum(WT_Downsampled, timep, freq, ...
                                                parameters.oddballTask.ERP_baselineCorrection, parameters.timeFreq.timeResolutionDivider, parameters, handles);
                            power = analyze_normalizeWaveletSpectrum(power, timep, freq, ...
                                                parameters.oddballTask.ERP_baselineCorrection, parameters.timeFreq.timeResolutionDivider, parameters, handles);
                        end

                    % Check if the interpolation failed
                    if noOfNaNs_afterInterpolation > noOfNaNs_beforeInterpolation
                        warning('Number of NaNs increased after interpolation!')
                    elseif noOfNaNs_afterInterpolation ~= 0
                        warning('Some NaNs after interpolation!')
                    else
                        WTout(ep,:,:) = WT_Downsampled;
                        powerOut(ep,:,:) = power;
                    end
                    
                    % EXTRA OUTPUTS (PARAMETERS CALCULATED)
                    [ITPC, ITLC, ITLCN, ERSP, avWT, WTav, avWTi, WTavi] = ...
                        analyze_wavelet_extraParametersAccum(WT_Downsampled, ITPC, ITLC, ITLCN, ERSP, avWT, WTav, avWTi, WTavi, parameters, handles);

                    % write out the results for better human readability of the code                   
                    % [noOfEpochs2, noOfScales, noOfTimePoints] = size(power);
                    averPowerPerChannel = squeeze(nanmean(powerOut(:,:,:),1));
                    stdOfPowerPerChannel = squeeze(nanstd(powerOut(:,:,:),1));
                    % [noOfScales, noOfTimePoints] = size(averPowerPerChannel);
                    % disp([noOfScales noOfTimePoint])

                end
                
                [extraParam.ITPC{ch}, extraParam.ITLC{ch}, extraParam.ITLCN{ch}, extraParam.ERSP{ch}, extraParam.avWT{ch}, extraParam.WTav{ch}, extraParam.Induced{ch}] = ...
                    analyze_wavelet_extraParameters(ITPC, ITLC, ITLCN, ERSP, avWT, WTav, avWTi, WTavi, noOfEpochs, parameters, handles);

                % now the different channels might have different amount of
                % epochs so we cannot put all the values to a single matrix, so
                % we use a cell
                %averPowerPerChannel = nanmean(WTout(:,:,:),1);
                %stdOfPowerPerChannel = nanstd(WTout(:,:,:),1);
                averageOfTheChannel{ch}(:,:) = averPowerPerChannel;
                stdOfTheChannel{ch}(:,:) = stdOfPowerPerChannel;

            end       

            fprintf('\n')
            disp(['            - timeDiv: ', num2str(parameters.timeFreq.timeResolutionDivider), ', min f: ', num2str(min(freq)), ', max f: ', num2str(max(freq)), ', freqRes: ', num2str(freq(2)-freq(1)), ' Hz'])       

    %% PLOT THE TIME-FREQUENCY 
    
        isNaN = logical(isNaN);        
        [noOfScales, noOfTimePoints, noOfEpochs] = size(WTout);                   
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
                [c, handle_contours{ch}] = contourf(T, F, Z, contourLevels, 'EdgeColor', 'none');                    
                epochLimits(3,ch,:) = [min(min(Z)) max(max(Z))];
                
                lab(ch,1,3) = xlabel('Time [ms]'); 
                lab(ch,2,3) = ylabel('Frequency [Hz]');                    
                tit(ch,3) = title(['n = ', num2str(sum(isNaN(:,ch) == 0)), ', IAF = ', num2str(IAF_peak), ' Hz']);
                ylim([1 30])

                % annotate the individual alpha range                    
                alphaLimits = IAF_peak + parameters.powerAnalysis.eegBins.freqs{1}; % 1 is now hard-coded, see init_defaultParameters.m to make this more robust

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
        set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold') 
        set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2, 'FontWeight', 'bold')            
        
            for ch = 1 : noOfChannels
                caxis(sp(3,ch), [min(min(epochLimits(3,:,:))) max(max(epochLimits(3,:,:)))]);
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
        


