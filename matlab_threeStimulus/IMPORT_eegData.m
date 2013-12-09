function [dataMatrix, triggers, info, sampleRate] = IMPORT_eegData(i, fileNameIn, inputFiles, handles)   

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempIMPORT.mat';
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

    % Raw import, import as it is saved
    % fileNameIn
    disp(['  Import the BDF file (', fileNameIn{i}, ')'])
    [info, dataMatrix, triggersRaw] = import_AndyReadBDFrev3(inputFiles, handles.parameters.SignalsToRead, handles.parameters.EEG.invertPolarity);     
    sampleRate = info.srate; % [Hz]
    info.calFactors = (info.PhysMax-info.PhysMin)./(info.DigMax-info.DigMin);
    
        if handles.flags.showDebugMessages == 1; 
            disp(info)            
        end
        
        if handles.parameters.EEG.invertPolarity == 1
           disp('INVERT POLARITY was TRUE so EEG signal is inverted!'); disp(' ');
        end

    % Process the triggers
    triggers = import_processTriggers(triggersRaw, handles.parameters.triggerSignals, handles.parameters.triggerPrecision, handles);
    
    % Now the EEG data file is started manually some seconds before the
    % PsychoPy routine is started, thus there is some extra at the end and
    % at the beginning of the .BDF file so we can trim this excess out
    % using the trigger .recON which is set high by the PsychoPy routine
    % when the experiment is started
    
        % first the data matrix (EEG + EOG + ECG)
        dataMatrixTemp = dataMatrix;
        dataMatrix = dataMatrixTemp(triggers.recON,:);
        
        % correct the triggers also
        triggerNames = fieldnames(triggers);
        
            for i = 1 : length(triggerNames)            
                triggers.(triggerNames{i}) = triggers.(triggerNames{i})(triggers.recON);                
            end
    
        samplesTrimmed = length(dataMatrixTemp) - length(dataMatrix);
        samplesTrimSec = samplesTrimmed / sampleRate;
        lengthSec = length(dataMatrix) / sampleRate;
        disp(['     .. trimmed off ', num2str(samplesTrimmed),  ' samples (', num2str(samplesTrimSec, 4), ' seconds) [.recON trigger]'])
        disp(['       .. duration of the recording: ', num2str(lengthSec, 4), ' seconds (', num2str(length(dataMatrix)), ' samples @ ', num2str(sampleRate), ' Hz)'])
    
    
