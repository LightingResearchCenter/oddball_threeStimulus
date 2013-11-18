function [epochs, stimulusON] = pre_findERP_Epochs(data, triggers, triggerUsed, stimulusType, alpha, found_oddballON, epochIndices_IN, baselineCorr, endOfEpochCorr, epochs, parameters, handles)

    
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if 1 == 1
        debugMatFileName = 'tempFindEpochs_ERPs.mat';
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

    j = 1;
        initCount = 0;
    irr = 1;
    reg = 1;
    
    % whos
    % length(data)
    % length(triggers.oddTone)    
    % CORRECT AT SOME POINT, GIVE THE TRIGGER TYPE as INPUT    

    debugEpochExtraction = 0;
    
    if strcmp(stimulusType, 'standard')
        divider = 40;
    else
        divider = 10;
    end
        
        
    if isempty(epochIndices_IN)
        fprintf(['      .. ', stimulusType,' (finding indices, slow): '])  
        % now we do the computationally heavy index finding
        
        for i = 2 : (length(data) - 1)            

            if initCount == 0 || initCount ~= j
                stimulusON(j,:) = [0 0]; % initialize
                initCount = j;                
            end

            if triggerUsed(i) == 1 && triggerUsed(i-1) == 0

                stimulusON(j,1) = i; % start     
                originalOnset(j) = stimulusON(j,1);

                    % Now the stimulusON found above corresponds to the trigger
                    % sent out by the Python script when we "wanted" the audio
                    % to be played, in reality, there is a delay to the actual
                    % onset of audio which can be found by finding the first
                    % onset after above-found index when the audio trigger is
                    % high (directly connected to the audio output)
                    searchWindowSecs = 1.5;
                    audioTriggerTemp = triggers.audio(stimulusON(j,1):stimulusON(j,1)+(parameters.EEG.srate*searchWindowSecs));
                    audioOnsetIndex = find(audioTriggerTemp == 0, 1, 'first');
                    audioOnsetIndex(j) = audioOnsetIndex + stimulusON(j,1); % add back the trimmed portion
                    delay(j) = audioOnsetIndex(j) - stimulusON(j,1);

                    % disp([stimulusON(j,1) audioOnsetIndex delay(j)]) 

                    %{
                    figure,                
                    hold on
                    plot(triggers.oddTone(stimulusON(j,1)-100:audioOnsetIndex+100), 'b')
                    plot(triggers.audio(stimulusON(j,1)-100:audioOnsetIndex+100), 'r')
                    %}
                                

                % correct stimulusOnset
                stimulusON(j,1) = audioOnsetIndex(j);
                    
                % add the preOnset baseline needed for ep_den
                stimulusON(j,1) = stimulusON(j,1) - baselineCorr;

                if strcmp(stimulusType, 'standard')

                    % now we have 240 standard tones compared to 40 odd/target
                    % and 40 distracter/non-targets

                    % found_oddballON

                end

            elseif triggerUsed(i) == 1 && triggerUsed(i+1) == 0 && stimulusON(j,1) ~= 0

                stimulusON(j,2) = i; % end            

                durationOfStimulusWithJitter(j) = stimulusON(j,2) - (originalOnset(j) + baselineCorr); % duration of oddball in samples            

                    % the oddball OFFSET is now found at 200 ms (as that is how
                    % the trigger is defined), so we need to add 0.6 sec *
                    % sampleRate to the oddball OFFSET
                        %toBeAdded = parameters.oddballTask.SOA_duration - parameters.oddballTask.triggerDuration;
                        %stimulusON(j,2) = stimulusON(j,2) + (toBeAdded * parameters.EEG.srate);

                    % you can later correct this for the actual audio ON or
                    % just manually add the required time for the stimulusON and
                    % assume that everything is properly phase-locked with the
                    % AUDIO ONSET
                        % stimulusON(j,2) = stimulusON(j,1) + ((parameters.oddballTask.ERP_duration + parameters.oddballTask.ERP_baseline) * parameters.EEG.srate);                                        

                % correct the end to be exactly the same for all the epoch as
                % there is small jitter in the trigger durations, redundancy
                % here as now number of samples just added to start point, but
                % this part of the code kept here to quantify the actual jitter
                % if needed/wanted
                stimulusON(j,2) = stimulusON(j,1) + endOfEpochCorr -1;                               
                    
                if rem(j,divider) == 0
                    fprintf('%s%s', num2str(j), ' ')
                end         
                
                % DEBUG
                %{
                a1 = stimulusON(j,1)
                a2 = stimulusON(j,2)
                size(data)
                plot(data(stimulusON(j,1):stimulusON(j,2),1:4))
                title(num2str(j))
                pause(0.3)
                %}

                j = j + 1; % increment the accumulator index
                
            end
        end
        
        % remove the last row that might be zeroes
        if stimulusON(end,1) == 0
            stimTemp = stimulusON;
            stimulusON = stimTemp(1:end-1,:);
        end
        
        % quantify delay between trigger and audio onset
        epochs.meanDelay = nanmean(delay);
        epochs.meanDelayStd = nanstd(delay);
        epochs.meanDelayMilliSec = 1000*epochs.meanDelay/parameters.EEG.srate;
        epochs.meanDelayMilliSecStd = 1000*epochs.meanDelayStd/parameters.EEG.srate;

        % quantify the jitter in epoch duration as defined from trigger onset
        % and offset
        epochs.meanStimulusDuration = nanmean(durationOfStimulusWithJitter);
        epochs.meanStimulusDurationStd = nanstd(durationOfStimulusWithJitter);
        epochs.meanStimulusDurationMilliSec = 1000*epochs.meanStimulusDuration/parameters.EEG.srate;
        epochs.meanStimulusDurationMilliSecStd = 1000*epochs.meanStimulusDurationStd/parameters.EEG.srate;
        
        fprintf('%s\n', ' ') % line change

        

        
    else
        % now we have done it already so we could just avoid the for-loop
        if strcmp(stimulusType, 'oddball')
            stimulusON = epochIndices_IN.oddballON;
        elseif strcmp(stimulusType, 'distracter')
            stimulusON = epochIndices_IN.distracterON;
        elseif strcmp(stimulusType, 'standard')
            stimulusON = epochIndices_IN.stdON;
        else
            error('Problem with the stimulusType variable, typo?')
        end
    end
    
        
    %% Now go through the found indices        
    fprintf(['        .. finding epochs, fast: '])  
    for j = 1 : length(stimulusON)
        
        indAudioOn = stimulusON(j,1) + baselineCorr;
        
        durationOfStimulus(j) = stimulusON(j,2) - stimulusON(j,1);
        epochsRAW.ERP{j} = data(stimulusON(j,1):stimulusON(j,2), :); % time domain EEG of the oddball        

        % if whole epoch is artifacted
        if sum(sum(isnan(epochsRAW.ERP{j}))) > 0
            epochsRAW.ERP{j}(:,:) = NaN;
            RMS_alphaFixed = NaN;
            RMS_alphaIAF = NaN;
            epochsRaw.RT(j) = Inf; % NaN if epoch is okay but subject failed to respond

        else           

            % remove the DC offset / drift / trend from the epoch
            epochsRAW.ERP{j} = pre_removeBaseline_epochByEpoch(epochsRAW.ERP{j}, j, parameters, handles);

            % Do power analysis of the baseline period before the
            % stimulus as done by Barry et al. (2000)
            [RMS_alphaFixed, RMS_alphaIAF] = pre_baselinePowerAnalysis(epochsRAW.ERP{j}, baselineCorr, alpha, parameters, handles);

            % vectors of time and button presses, if you need to output,
            % index these
            reactionRaw.timeVec = (-durationOfStimulus(j):1:durationOfStimulus(j))';
            reactionRaw.button  = triggers.button(stimulusON(j,1) : (stimulusON(j,2)+durationOfStimulus(j)));

            % DEBUG        
            if (j == 241 || (length(stimulusON) == 39 && j == 39) || j == 41) && debugEpochExtraction == 1

                offSetPadding = 128;

                ind1 = stimulusON(j,1) - offSetPadding;
                ind2 = stimulusON(j,2) + offSetPadding;

                %length(data)
                %whos
                if ind2 > length(data)
                    ind2 = length(data)
                end
                
                dataIn = data(ind1:ind2,:);

                button_debug = triggers.button(ind1:ind2);
                audio_debug = triggers.audio(ind1:ind2);
                trigger_debug = triggerUsed(ind1:ind2);

                x_debug = (1 : 1 : length(button_debug))' / parameters.EEG.srate;
                x_debug = linspace(stimulusON(j,1), stimulusON(j,2), length(button_debug))';

                whos            

                subplot(2,1,1)
                    plot(x_debug, dataIn)

                subplot(2,1,2)
                plot(x_debug, button_debug, 'r', x_debug, audio_debug, 'b', x_debug, trigger_debug, 'k')
                    legend('Button', 'Audio', 'Python Trigger',3)
                    legend('boxoff')


                title(num2str(j))
                drawnow()
                pause(1.0)            

            end        

            % Get the reaction time in seconds, note that the reaction time
            % could be negative as the SOA is fixed and most likely some
            % people start anticipating the following std tone
            indButton = stimulusON(j,1) + find(reactionRaw.button == 1, 1, 'first');
            reactionTimeInSamples = indButton - indAudioOn;

            if ~isempty(indButton)
                epochsRaw.RT(j) = reactionTimeInSamples / parameters.EEG.srate;
                if epochsRaw.RT(j) > parameters.oddballTask.SOA_duration
                    epochsRaw.RT(j) = NaN;
                    % make sure that no button presses are detected after that
                    % epoch
                end
            else
                epochsRaw.RT(j) = NaN; % lapse, no reaction time found during the following std tone
            end
        
        
        %% ADD SOME FIELDS TO RT later
            % .latency
            % .duration

        %% ADD DEBUG 
            % .triggerDuration
            % .trigger-audio delay
            % .audioDuration
            
        end
        %}

        
        % there used to be a switch between irregular and regular
        % target-to-target interval so that is why we kinda re-assign
        % the same variable to something else
        epochs.ERP{j} = epochsRAW.ERP{j};
        epochs.RT(j) = epochsRaw.RT(j);
        epochs.alphaFixed(j) = RMS_alphaFixed;
        epochs.alphaIAF(j) = RMS_alphaIAF;                    

        if rem(j,divider) == 0
            fprintf('%s%s', num2str(j), ' ')
        end
         
        
    end
    fprintf('%s\n', ' ') % line change    
    
    if debugEpochExtraction == 1
        % debug_plotEPochExctraction(epochs.ERP, epochs.alphaIAF, baselineCorr, endOfEpochCorr, stimulusType, parameters, handles)
    end
    
    