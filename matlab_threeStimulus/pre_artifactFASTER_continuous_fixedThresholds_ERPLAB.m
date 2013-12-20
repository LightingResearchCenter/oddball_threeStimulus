function [NaN_indices_EEG, NaN_indices_EOG, indices_moving, indices_movingEOG, indices_step, vDiffOutMovWindow, vDiffOutMovWindowEOG, vDiffOutStep, isNaN] = ...
    pre_artifactFASTER_continuous_fixedThresholds_ERPLAB(EEG, EOG, ECG, debugOn, erpType, parameters, handles)

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempFASTERfixedThresholdsContinuous.mat';
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
        
    [numberOfChannels, noOfSamples, numberOfEpochs] = size(EEG);
    NaN_indices_EEG = zeros(numberOfEpochs, parameters.EEG.nrOfChannels);
    NaN_indices_EOG = zeros(numberOfEpochs, parameters.EEG.nrOfChannels);
    
    for ep = 1 : numberOfEpochs                   

                    
        %% ERPLAB ARTIFACT DETECTION          
                   
           [indices_moving(ep,:), indices_movingEOG(ep,:), indices_step(ep,:), vDiffOutMovWindow(ep,:,:), vDiffOutMovWindowEOG(ep,:,:), vDiffOutStep(ep,:,:)] = ...
                pre_artifact_ERPLAB(squeeze(EEG(:,:,ep)), squeeze(EOG(:,:,ep)), ep, parameters);            

            
        %% whether EEG channels exceeds the threshold
            
            % "Stupid" fixed threshold rejection  
            maxAmplitudes = max(abs(EEG(:,:,ep)));
            NaN_indices_EEG_raw = abs(EEG(:,:,ep)) > parameters.artifacts.fixedThr; % individual samples                        
            NaN_indices_EEG(ep,:) = logical(sum(NaN_indices_EEG_raw,2)); % for channel
            
            %disp(maxAmplitudes)
            %disp(NaN_indices_EEG(ep,:))            
            
            % whether EOG channel exceeds the threshold
            NaN_indices_EOG_raw = abs(EOG(:,:,ep)) > parameters.artifacts.fixedThrEOG; % individual samples
            NaN_indices_EOG_raw = logical(sum(NaN_indices_EOG_raw)); % for EOG channel        
            NaN_indices_EOG(ep,:) = repmat(NaN_indices_EOG_raw, length(NaN_indices_EOG_raw), length(NaN_indices_EEG(ep,:))); % to match the EEG channels
         
        
        
    end
    
    % replicate vectors to all channels
    %indices_moving = repmat(indices_moving', 1, parameters.EEG.nrOfChannels);
    indices_step = repmat(indices_step, 1, parameters.EEG.nrOfChannels);
    
    % simple boolean, giving 1 when the epoch is artifacted 
    isNaN = logical(NaN_indices_EEG + NaN_indices_EOG + indices_moving + indices_step);

    %% DEBUG PLOTS

        %% EOG TIME DOMAIN PLOT
        if debugOn == 2

            scrsz = get(0,'ScreenSize'); % get screen size for plotting
            fig = figure('Color', 'w', 'Name', ['Continuous FASTER']);
                set(fig, 'Position', [0.35*scrsz(3) 0.04*scrsz(4) 0.45*scrsz(3) 0.90*scrsz(4)])
                
            
            EEGmat = pre_epochMatrixToVector(EEG, parameters, handles);
            EOGmat = pre_epochMatrixToVector(EOG, parameters, handles);
            rows = 4; cols = 1;

            ind = 1;
            sp(ind) = subplot(rows,cols,ind);
                p(ind,:) = plot(EEGmat)

            ind = 2;
            sp(ind) = subplot(rows,cols,ind);
                p(ind,:) = plot(EOGmat)

        end
        
        
    

%% SUBFUNCTION WRAPPER FOR THE CRAP from ERPLAB
function [moving_isNaN, movingEOG_isNaN, step_isNaN, vDiffOutMovWindow, vDiffOutMovWindowEOG, vDiffOutStep] = pre_artifact_ERPLAB(epochIn, EOG, ep, parameters)

    

    % Use EEGLAB structure fields
    EEG.data = epochIn;
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

        EEG.data = EOG;
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
        EEG.data = EOG;
        

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