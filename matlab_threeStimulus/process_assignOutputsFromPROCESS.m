function process_assignOutputsFromPROCESS(epochs, ERP_components, alpha, powers, SEM, heart, EOG, fractalAnalysis, timeFreqEpochs, timeFreq, TF_allEpochs, TF_derivedMeasures, parameters, handles)


    % you can discard some data, and make sure that you only the needed
    % data is saved to disk, or if you just want to optimize the use of
    % disk space or something    
    
    % Also it is nicer to save the EPOCHS and timeFreq data to its own .mat
    % file from ERP components as it is easier to batch plot them then
    % later as you don't have to load such big chunks of data at once
    % making processing faster (as well debugging and development
    
    
    %% DEBUG
    debugMatFileName = 'tempAssignOutput.mat';
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
    whos
    
    
    % get rid of the "double existence" of the parameters structures
    parameters = handles.parameters;
    handles.parameters = [];

    %% ERP Waveforms

        fileMatOut = strrep(handles.inputFile, '.bdf', '_fullEpochs.mat');        
        disp(['      .. All ERP Epochs [', fullfile(handles.path.matFilesOut, fileMatOut), ']'])
        save(fullfile(handles.path.matFilesOut, fileMatOut), 'epochs', 'parameters', 'handles')        
        
    %% ERP Waveforms stat, instead of outputting all the epochs, calculate the mean/median/etc. here and save
    
        fileMatOut = strrep(handles.inputFile, '.bdf', '_statEpochs.mat');        
        disp(['      .. Stats of ERP Epochs [', fullfile(handles.path.matFilesOut, fileMatOut), ']'])
        
        % quick'n'dirty fix
        erpFilterType = {'bandpass'};
        erpType = {'ERP'};
        chsToPlot = {'Cz'; 'Pz'};
        responseTypes = fieldnames(epochs.(erpFilterType).(erpType));
        
        for filt = 1 : length(erpFilterType)       
            for erp = 1 : length(erpType)  
                for resp = 1 : length(responseTypes)            
                    % process using a subfunction
                    [epochsStat, debug] = batch_timeDomEpochsAverage(epochs.(erpFilterType).(erpType).(responseTypes{resp}), thisFileMarkedAsOutlier, handles.parameters, handles);
                end
            end
        end
        
        % save
        save(fullfile(handles.path.matFilesOut, fileMatOut), 'epochsStat', 'parameters', 'handles')        

    %% Derived ERP Measures such as amplitude latency etc.

        % "Auxiliary measures" (scalars)            
        analyzed_aux.IAF_amplitGravity = alpha.amplit_gravity;
        for i = 1 : length(parameters.powerAnalysis.eegBins.freqs)
            fieldName = powers{i}.label;
            analyzed_aux.Amplit.(fieldName) = powers{i}.ps.powerData;                
            analyzed_aux.PSD.(fieldName) = powers{i}.PSD.powerData;
        end

        fileMatAnalyzedOut = strrep(handles.inputFile, '.bdf', '_analyzed.mat');                        
        disp(['       .. Auxiliary measures [', fullfile(handles.path.matFilesOut, fileMatOut), ']'])
        save(fullfile(handles.path.matFilesOut, fileMatAnalyzedOut), 'ERP_components', 'analyzed_aux', 'powers', 'alpha', 'parameters', 'handles') % need powers and alpha?

    %% Time-Frequency derived
        fileMatOutWavelet = strrep(handles.inputFile, '.bdf', '_waveletDerived.mat');
        disp(['        .. Derived Time-Frequency measures [', fullfile(handles.path.matFilesOut, fileMatOutWavelet), ']'])
        save(fullfile(handles.path.matFilesOut, fileMatOutWavelet), 'timeFreqEpochs', 'timeFreq', 'TF_derivedMeasures', 'parameters', 'handles')


    %% Time-Frequency, all epochs 
        % (REQUIRE A LOT OF DISK SPACE), around ~1 GB of session -> ~120 GB for 10 subjects
        %{
        fileMatOutWaveletEpochs = strrep(handles.inputFile, '.bdf', '_waveletEpochs.mat');
        disp(['       .. All Wavelet (TF) Epochs) [', fullfile(handles.path.matFilesOut, fileMatOutWaveletEpochs), ']'])
        save(fullfile(handles.path.matFilesOut, fileMatOutWaveletEpochs), 'TF_allEpochs', 'handles')
        %}

    %% EEG Fractal Analysis
        fileMatOutFractalEEG = strrep(handles.inputFile, '.bdf', '_fractalEEG.mat');
        disp(['         .. EEG Fractal Analysis [', fullfile(handles.path.matFilesOut, fileMatOutFractalEEG), ']'])
        save(fullfile(handles.path.matFilesOut, fileMatOutFractalEEG), 'fractalAnalysis', 'parameters', 'handles')

    %% "Extra measures", such as EOG and HEART data
        analyzed_extraSensors.SEM = SEM;
        analyzed_extraSensors.heart = heart;

        fileMatOutExtra = strrep(handles.inputFile, '.bdf', '_extraAnalysis.mat');
        disp(['          .. Extra sensors (EOG,ECG) [', fullfile(handles.path.matFilesOut, fileMatOutExtra), ']'])
        save(fullfile(handles.path.matFilesOut, fileMatOutExtra), 'analyzed_extraSensors', 'parameters', 'handles')


