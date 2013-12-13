function [timeFreqEpochs, timeFreq_target, timeFreq_distr, timeFreq_std, allEpochs, derivedMeasures] = ...
                pre_waveletWrapperForAllStimuli(epochs_raw, artifactIndices_FASTER, epochs_distr_raw, artifactIndices_FASTER_distr, epochs_std_raw, artifactIndices_FASTER_std, ...
                IAF, parameters, handles)
            
    % Just so slow to run the PROCESS_singleFile from the start, so not the
    % most elegant wrapper necessarily, but just easier to debug
    
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempWaveletWrapper.mat';
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
    
    % TARGET
    [timeFreqEpochs.target.real, timeFreqEpochs.target.imag, allEpochs.target, derivedMeasures.target, timeFreq_target] = ...
        analyze_timeFreqAnalysisForCondition(epochs_raw, artifactIndices_FASTER, 'target', IAF, parameters, handles);   
    
    % DISTRACTER
    [timeFreqEpochs.distr.real, timeFreqEpochs.distr.imag, allEpochs.distr, derivedMeasures.distr, timeFreq_distr] = ...
        analyze_timeFreqAnalysisForCondition(epochs_distr_raw, artifactIndices_FASTER_distr, 'distracter', IAF, parameters, handles);
    
    % STANDARD
    [timeFreqEpochs.std.real, timeFreqEpochs.std.imag, allEpochs.std, derivedMeasures.std, timeFreq_std] = ...
        analyze_timeFreqAnalysisForCondition(epochs_std_raw, artifactIndices_FASTER_std, 'standard', IAF, parameters, handles);