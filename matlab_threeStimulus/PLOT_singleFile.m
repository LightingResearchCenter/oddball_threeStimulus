function PLOT_singleFile(epochs_target, epochs_distracter, epochs_standard, ...
             epochsEP_target, epochsEP_distracter, epochsEP_standard, ...
             analyzed_target, analyzed_distracter, analyzed_standard, ...
             analyzed_aux, analyzed_extraSensors, analyzed_fractal, analyzed_TF, ...
             dataMatrix_filtGeneral, alpha, powers, info, handles)

    debugMatFileName = 'tempPlotWaveforms.mat';
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
    
    % INPUT EPOCHS
    % epochs_Jongsma
    % epochs_Jongsma_CNV
    % epochs_Jongsma_filt
    % epochs_Jongsma_raw
    % epochs_raw          - 1) Raw EEG without bandpass filtering, artifact rejection
    % epochs_filt         - 2) Bandbass and artifact rejected epochs
    % epochs_ep           - 3) EP Denoised epochs after _filt 

    disp(['  Plotting the data'])
    numberOfSamplesPerEpoch = abs((-handles.parameters.oddballTask.ERP_duration - handles.parameters.oddballTask.ERP_baseline)) * handles.parameters.EEG.srate;
    handles.parameters.oddballTask.timeVector = linspace(-handles.parameters.oddballTask.ERP_baseline, handles.parameters.oddballTask.ERP_duration, numberOfSamplesPerEpoch)';

    %% WAVEFORM PLOT
    
        % plot_timeFrequency(analyzed_TF, handles)
        
        plot_waveforms(epochs_target, epochs_distracter, epochs_standard, ...
             epochsEP_target, epochsEP_distracter, epochsEP_standard, ...
             analyzed_target, analyzed_distracter, analyzed_standard, handles)                
    
        % plot_componentTrend(analyzed_target, analyzed_distracter, analyzed_standard, 'allSamples', handles)        
        
        % plot_componentTrend(analyzed_target, analyzed_distracter, analyzed_standard, 'averaged', handles)