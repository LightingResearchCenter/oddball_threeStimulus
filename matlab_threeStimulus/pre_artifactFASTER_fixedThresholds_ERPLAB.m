function [NaN_indices_EEG, NaN_indices_EOG, indices_moving, indices_movingEOG, indices_step, vDiffOutMovWindow, vDiffOutMovWindowEOG, vDiffOutStep, isNaN] = ...
    pre_artifactFASTER_fixedThresholds_ERPLAB(EEG, EOG, ECG, debugOn, erpType, parameters, handles)

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempFASTERfixedThresholds.mat';
        if nargin == 0
            load('debugPath.mat')
            load(fullfile(path.debugMATs, debugMatFileName))
            close all
        else
            if handles.flags.saveDebugMATs == 1
                if ~strcmp(erpType, 'fff') % 'standard')
                    % do not save for standard tone as there are so many
                    % trials that debugging and developing of this function
                    % is so much slower compared to target and distracter
                    path = handles.path;
                    save('debugPath.mat', 'path')
                    save(fullfile(path.debugMATs, debugMatFileName))            
                end
            end
        end 
    end
        
    numberOfEpochs = length(EEG.ERP);
    
    NaN_indices_EEG = zeros(numberOfEpochs, parameters.EEG.nrOfChannels);
    NaN_indices_EOG = zeros(numberOfEpochs, parameters.EEG.nrOfChannels);
    
    HDR.SampleRate = parameters.EEG.srate;
    
    % if you use the very raw input without referencing
    % parameters.BioSemi.chOffset = 2;

    if parameters.artifacts.FASTER_applyRegressionEOG == 1                    
        fprintf(' regression-EOG ');
    end

    if parameters.artifacts.FASTER_applyRegressionECG == 1
        fprintf(' regression-ECG (PRA) ');
    end
    
    % disp('        - Epoch: ')
    for ep = 1 : numberOfEpochs    
        
        %% EEG
        
            % slightly redundant processing (fix maybe later for computational
            % efficiency)
            noOfChannels = size(EEG.ERP{ep},2);
            if noOfChannels ~= 6
                Reference = mean(EEG.ERP{ep}(:,1:parameters.BioSemi.chOffset),2);
                EEG_channels = EEG.ERP{ep}(:,parameters.BioSemi.chOffset+1:parameters.BioSemi.chOffset+parameters.EEG.nrOfChannels);            
                    EEG_channels = EEG_channels - repmat(Reference, 1, size(EEG_channels,2));
            else
                EEG_channels = EEG.ERP{ep}(:,1:parameters.EEG.nrOfChannels);
            end

                noOfChannels = size(EEG_channels,2);
                for ch = 1 : noOfChannels
                    % EEG_channels(:,ch) = pre_remove5060hz_modified(EEG_channels(:,ch), HDR, 'PCA 60');                
                end
                
            % EOG
            EOGin = EOG.ERP{ep}(:,1);
                % EOG already referenced so no need to subtract

                % filter            
                EOGnotch = pre_remove5060hz_modified(EOGin, HDR, 'PCA 60');

                cutOffs = [45 0]; % remove high frequency noise, try to keep the trend
                filterOrder = 10; % conservative 
                try 
                    EOGnotch = pre_bandbassFilter(EOGnotch, HDR.SampleRate, cutOffs, filterOrder, filterOrder, handles);
                catch err
                    err
                    filterOrder = filterOrder / 2;
                    EOG = pre_bandbassFilter(EOGnotch, HDR.SampleRate, cutOffs, filterOrder, filterOrder, handles);
                    error('Bandpass filter parameters too "extreme"')
                end
                
            % ECG
            ECGin = ECG.ERP{ep}(:,1);
                % NOTE! ECG unfiltered, maybe filter later, as now not used

            % Regress_EOG for EOG and ECG if wanted
                if parameters.artifacts.FASTER_applyRegressionEOG == 1                    
                    EEG_channels = pre_artifactRegressionWrapper(EEG_channels, EOGnotch, parameters.EEG.nrOfChannels, parameters);                    
                end

                if parameters.artifacts.FASTER_applyRegressionECG == 1                    
                    EEG_channels = pre_artifactRegressionWrapper(EEG_channels, ECGin, parameters.EEG.nrOfChannels, parameters);
                end
                
             % re-correct for baseline        
            EEG_channels = pre_removeBaseline_epochByEpoch(EEG_channels, ep, parameters, handles);   
            

                    
        %% ERPLAB ARTIFACT DETECTION          
                   
           [indices_moving(ep,:), indices_movingEOG(ep,:), indices_step(ep,:), vDiffOutMovWindow(ep,:,:), vDiffOutMovWindowEOG(ep,:,:), vDiffOutStep(ep,:,:)] = ...
               pre_artifact_ERPLAB(EEG_channels, EOGnotch, ep, parameters);            
            
            % debugging            
            if ep == 40
               %{
               figure('Name', 'CRAP Debug')
               subplot(4,1,1)
               y = squeeze(vDiffOutMovWindow(ep,:,:))';
               plot(y); title(['MovingWindow (CRAP), max Diff = ', num2str(max(y))]);
               
               subplot(4,1,2)
               plot(EEG_channels); title('EEG');
               
               subplot(4,1,3)
               y = squeeze(vDiffOutMovWindowEOG(ep,:,:))';
               plot(y); title(['MovingWindow EOG (CRAP), max Diff = ', num2str(max(y))]);
               
               subplot(4,1,4)
               plot(EOGnotch); title('EOG');
               drawnow
               %fd
               %}
            end
            
            
        %% whether EEG channels exceeds the threshold
            
            % "Stupid" fixed threshold rejection  
            maxAmplitudes = max(abs(EEG_channels));
            NaN_indices_EEG_raw = abs(EEG_channels) > parameters.artifacts.fixedThr; % individual samples            
            NaN_indices_EEG(ep,:) = logical(sum(NaN_indices_EEG_raw,1)); % for channel
            
            %disp(maxAmplitudes)
            %disp(NaN_indices_EEG(ep,:))
            
            
            % whether EOG channel exceeds the threshold
            NaN_indices_EOG_raw = abs(EOGnotch) > parameters.artifacts.fixedThrEOG; % individual samples
            NaN_indices_EOG_raw = logical(sum(NaN_indices_EOG_raw)); % for EOG channel        
            NaN_indices_EOG(ep,:) = repmat(NaN_indices_EOG_raw, length(NaN_indices_EOG_raw), length(NaN_indices_EEG(ep,:))); % to match the EEG channels
         
        %% DEBUG PLOTS

            %% EOG TIME DOMAIN PLOT
            if debugOn == 1   
                %{
                try
                    subplot(5,8,ep)
                    hold on
                    plot(vDiffOutMovWindow(ep,:)); 
                    plot(vDiffOutStep(ep,:), 'k');
                    hold off
                    title(num2str(ep)); lab(1,ep) = xlabel('Epoch samples'); lab(1,ep) = ylabel('\muV');
                    xlim([1 length(EOG)])
                    set(gca, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
                    set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
                catch err
                    err
                end
                %}
                
            end
            
        
            %% EOG TIME DOMAIN PLOT
            if debugOn == 1   
                %{
                try
                    subplot(5,8,ep)
                    plot(EOG); title(num2str(ep)); lab(1,ep) = xlabel('Epoch samples'); lab(1,ep) = ylabel('\muV');
                    xlim([1 length(EOG)])
                    set(gca, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
                    set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
                catch err
                    err
                end
                %}
            end
            
            %% EEG TIME DOMAIN PLOT
            if debugOn == 1     
                %{
                try
                    subplot(5,8,ep)
                    plot(EEG_channels); title(num2str(ep)); lab(1,ep) = xlabel('Epoch samples'); lab(1,ep) = ylabel('\muV');
                    xlim([1 length(EEG_channels)])
                    set(gca, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
                    set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
                catch err
                    err
                end
                drawnow
                %{
                if ep == numberOfEpochs
                    leg = legend('Fz', 'Cz', 'Oz', 'Pz');
                    set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
                    legend('boxoff')
                end
                %}
                %}
            end            
            
        
            %% FREQUENCY ANALYSIS of EPOCH       
            if debugOn == 1
                %{
                subplot(1,2,1)
                    plot(EOGin)
                    title(num2str(ep)); lab(1,1) = xlabel('Epoch samples'); lab(1,2) = ylabel('\muV');
                    xlim([1 length(EOG)])
                    set(gca, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
                    set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
                subplot(1,2,2)
                    % FFT
                    Fs = parameters.EEG.srate;   
                    freqRange = 0 : 0.1 : parameters.filter.bandPass_hiFreq;
                    nOverlap = parameters.powerAnalysis.nOverlap; % in %  
                    [f, EOGStruct.amplit, EOGStruct.PSD, EOGStruct.amplit_matrix, EOGStruct.PSD_matrix] = ...
                        analyze_powerSpectrum(EOGin, Fs, 'Tukey', 0.10, length(EOGin), length(EOGin), freqRange, nOverlap, 'EOG');
                    plot(f, 10*log10(EOGStruct.PSD.mean), 'b'); lab(2,1) = xlabel('Freq'); lab(2,2) = ylabel('Amplitude');
                    xlim([0.1 2*parameters.filter.bandPass_hiFreq])

                    %HDR.SampleRate = parameters.EEG.srate;
                    %EOG = pre_remove5060hz_modified(EOG, HDR, 'PCA 60');
                    [f, EOGStruct.amplit, EOGStruct.PSD, EOGStruct.amplit_matrix, EOGStruct.PSD_matrix] = ...
                        analyze_powerSpectrum(EOG, Fs, 'Tukey', 0.10, length(EOG), length(EOG), freqRange, nOverlap, 'EOG');
                    hold on
                    plot(f, 10*log10(EOGStruct.PSD.mean), 'r');
                    hold off
                    legend('EOG in', 'EOG with NOTCH60')

                    subplot(1,2,1)
                    hold on
                    plot(EOG, 'r')
                    hold off

                pause(1.0)
                %}
            end
        
    end
    
    % replicate vectors to all channels
    %indices_moving = repmat(indices_moving', 1, parameters.EEG.nrOfChannels);
    indices_step = repmat(indices_step, 1, parameters.EEG.nrOfChannels);
    
    % simple boolean, giving 1 when the epoch is artifacted 
    isNaN = logical(NaN_indices_EEG + NaN_indices_EOG + indices_moving + indices_step);

    

%% SUBFUNCTION WRAPPER FOR THE CRAP from ERPLAB
function [moving_isNaN, movingEOG_isNaN, step_isNaN, vDiffOutMovWindow, vDiffOutMovWindowEOG, vDiffOutStep] = pre_artifact_ERPLAB(epochIn, EOG, ep, parameters)

    % Use EEGLAB structure fields
    EEG.data = epochIn';
    EEG.epoch = [];
    EEG.srate = parameters.EEG.srate;
    EEG.pnts = length(epochIn);
    EEG.event(1).type = 0;
    EEG.xmin = 0;
    EEG.xmax = length(EEG.data) / parameters.EEG.srate;
    EEG.setname = 'dummy blank string';
    
    %{
    plot(epochIn)
    hold on
    plot(EOG, 'm')
    hold off
    pause(1.0)
    %}

    % See http://erpinfo.org/erplab/erplab-documentation/manual_4/Artifact_Detection.html

    %% The Moving Window Peak-to-Peak Function

        % The most broadly useful artifact detection function supplied by ERPLAB 
        % is the moving window peak-to-peak threshold function.  
        % Peak-to-peak amplitude is the difference between the most positive and 
        % most negative voltages within a window.  A moving window peak-to-peak amplitude 
        % function computes the peak-to-peak amplitude within a series of windows within each epoch. 
        [indices_moving, vDiffOutMovWindow] = crap_mod(EEG, parameters.artifacts.CRAP.movWind_ampTh, parameters.artifacts.CRAP.movWind_windowWidth, parameters.artifacts.CRAP.movWind_windowStep, ...
                                   1:parameters.EEG.nrOfChannels);
             
        moving_isNaN = logical(sum(indices_moving));
        

    %% Moving window for EOG channel

        EEG.data = EOG';

        % The step function was designed to find the step-like changes in voltage 
        % that are produced when subjects make saccadic eye movements, 
        % but it is useful for detecting other kinds of artifacts as well (e.g., blinks).  
        % The step function begins by defining a step-shaped function of a particular width 
        % (e.g., 200 ms at one voltage and then 200 ms at a different voltage)
        [indices_movingEOG, vDiffOutMovWindowEOG] = crap_mod(EEG, parameters.artifacts.CRAP.movWindEOG_ampTh, parameters.artifacts.CRAP.movWindEOG_windowWidth, parameters.artifacts.CRAP.movWindEOG_windowStep, ...
                                                1:1);

        movingEOG_isNaN = logical(sum(indices_movingEOG));
        
    %% STEP?
    
        % check if really correct
        EEG.data = EOG';

        % Probably use some other function than the crap_mod, check later?
        % http://erpinfo.org/erplab/erplab-documentation/manual_4/Artifact_Detection.html


        [indices_step, vDiffOutStep] = crap_mod(EEG, parameters.artifacts.CRAP.step_ampTh, parameters.artifacts.CRAP.step_windowWidth, parameters.artifacts.CRAP.step_windowStep, ...
                                                1:1);

        step_isNaN = logical(sum(indices_step));
        

        %{
        handles.parameters.artifacts.CRAP.step_ampTh = 10;
        handles.parameters.artifacts.CRAP.step_windowWidth = 50;
        handles.parameters.artifacts.CRAP.step_windowStep = 15;
        %}


    %% Debug plot

    %{
        ind = 1;
        sp(ind) = subplot(3,1,ind);                
            cla
            plot(epochIn)
            hold on
            plot(EOG, 'k')
            hold off
            title(['Epoch #', num2str(ep)])

        ind = 2;
        sp(ind) = subplot(3,1,ind);                
            cla
            plot(vDiffOutMovWindow)
            hold on
            line([1 length(epochIn)], [max(parameters.artifacts.CRAP.movWind_ampTh) max(parameters.artifacts.CRAP.movWind_ampTh)])
            title(['Moving Window'])

        ind = 3;
        sp(ind) = subplot(3,1,ind);                
            cla
            plot(vDiffOutStep)
            hold on
            line([1 length(epochIn)], [parameters.artifacts.CRAP.step_windowStep parameters.artifacts.CRAP.step_windowStep])
            title(['Step'])

        set(sp, 'XLim', [1 length(epochIn)])

            pause
    %}