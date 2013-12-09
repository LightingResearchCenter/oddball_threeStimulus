function [epochs, epochs_distr, epochs_std, epoch_indicesOut] = pre_epochToERPs(data, triggers, epochIndices_IN, alpha, parameters, dataType, handles)

    %{
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempEpochToERPs.mat';
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
    end
    %}
    
    
    % whos
    baselineCorr = (parameters.oddballTask.ERP_baseline * parameters.EEG.srate);
    endOfEpochCorr = baselineCorr + (parameters.oddballTask.ERP_duration * parameters.EEG.srate);
    
        if parameters.oddballTask.ERP_baseline ~= parameters.oddballTask.ERP_duration
            %warning('Your baseline and post-stimulus ERP windows are not the same length, considering making them equal for optimal performance of the EP_DEN denoising (init_defaultParameters)')
            % The EP Wavelet Denoising algorithm would require the baseline
            % and ERP duration to have the same length
            
            % Now actually best probably to fix the baselines just to be
            % different for EP DENOISING than for the rest of the analysis
            % script, time-frequency analysis requires rather long duration
            % after stimulus onset (900 ms in Peng et al. (2012), 
         
        end
    
        integerTest = mod(baselineCorr,1);
        if integerTest ~= 0
            %warning('non-integer epoch definition, check .ERP_baseline from init_DefaultParameters.m. Floored the index')
            baselineCorr = floor(baselineCorr);
        end

        integerTest = mod(endOfEpochCorr,1);
        if integerTest ~= 0
            %warning('non-integer epoch definition, check .ERP_duration from init_DefaultParameters.m. Floored the index')
            endOfEpochCorr = floor(endOfEpochCorr);
        end
    
    epochs.samplesPerEpoch = endOfEpochCorr;
    epochs_std.samplesPerEpoch = endOfEpochCorr;
    
    
    %% ODDBALLS           
    stimulusType = 'oddball';
    [epochs, oddballON] = pre_findERP_Epochs(data, triggers, triggers.oddTone, stimulusType, alpha, [], epochIndices_IN, baselineCorr, endOfEpochCorr, epochs, parameters, handles);
    
    %% DISTRACTERS
    stimulusType = 'distracter';
    [epochs_distr, distracterON] = pre_findERP_Epochs(data, triggers, triggers.distracter, stimulusType, alpha, [], epochIndices_IN, baselineCorr, endOfEpochCorr, epochs, parameters, handles);
            
    %% STANDARD TONES    
    stimulusType = 'standard';
    [epochs_std, stdON] = pre_findERP_Epochs(data, triggers, triggers.stdTone, stimulusType, alpha, oddballON, epochIndices_IN, baselineCorr, endOfEpochCorr, epochs, parameters, handles);
            
    
    if ~isempty(epochIndices_IN)
        epoch_indicesOut = epochIndices_IN;
    else % for the first run
        epoch_indicesOut.oddballON = oddballON;
        epoch_indicesOut.distracterON = distracterON;
        epoch_indicesOut.stdON = stdON;
    end
    
    
    % epochsRAW not returned now
    if strcmp(dataType, 'filt') == 1 % show only once
        try             
            reactionTimes = [epochs.RT'];
            
        catch err

            disp(err.identifier)
            warning('Dimensions are not the same for RTs of irregular and regular condition, you may want to check out why is this so?')
            RTs = epochs.RT';
            % irregRTs = epochs.RT_irregular';
            disp(['        Length of regular RTs: ', num2str(length(RTs))])
            % disp(['        Length of irregular RTs: ', num2str(length(irregRTs))])

            % This might occur by mistake if you have the EEG recording
            % ON and you initiate the task from PsychoPy and stop it
            % and do not open a new file for example
            % [oddballON irregularOn]
            %                 ans =
            % 
            %       240381      244476           1
            %       314243      318338           1
            %       350255      354350           1
            %       366622      370717           1
            %       392812      396907           1
            %       425547      429642           1
            %       458286      462381           1
            %       477926      482021           1
            %       500843      504938           1
            %       530653      534748           0

            % You might want to add something systematic here that
            % either checks the number of consecutive irregularOn
            % indices (9 in the above case instead of 8) or check the
            % interval between thee oddballs (the first one above
            % occurs a lot earlier than the "actual" oddballs in
            % relation to each other)
            disp('         .. Manual quick and dirty fix now, first value eliminated')
                for i = 2 : length(epochs.oddball_irregular)
                    epochs.oddball_irregularTemp{i-1} = epochs.oddball_irregular{i};                
                    epochs.RT(i-1) = epochs.RT(i);
                end
                epochs.oddball = epochs.oddball_irregularTemp;
                epochs.RT_irregular = epochs.RT_irregularTemp;

        end            

        % DEBUG INFO
        if handles.flags.showDebugMessages == 1
            disp('From: pre_epochToERPs.m')
                % durationOfOddball
                % durationOfOddballWithJitter
                % [oddballON irregularOn]          
                % reactionTimes = [epochs.RT_regular' epochs.RT_irregular']
                % whos
        end
    end
    
    if strcmp(dataType, 'filt') == 1 % show only once
        disp(' ')
        disp(['          .. Mean oddball trigger duration: ', num2str(nanmean(epochs.meanStimulusDuration)), ' +/- ', num2str(epochs.meanStimulusDurationStd), ' samples'])
        disp(['           .. Mean oddball trigger duration: ', num2str(epochs.meanStimulusDurationMilliSec), ' +/- ', num2str(epochs.meanStimulusDurationMilliSecStd), ' milliseconds, ideal theoretical = ', num2str(parameters.oddballTask.triggerDuration*1000), ' ms, NEGATIVE?'])
        disp(['            .. Mean delay (audio onset - trigger onset): ', num2str(nanmean(epochs.meanDelay)), ' +/- ', num2str(epochs.meanDelayStd), ' samples'])
        disp(['             .. Mean delay: ', num2str(epochs.meanDelayMilliSec), ' +/- ', num2str(epochs.meanDelayMilliSecStd), ' milliseconds'])
    end
    
    