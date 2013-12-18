function batchCompare_epochWaveforms(fileNameFields, outlierFilenameList, epochType, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempEpochCompare.mat';
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
    
    % define what you want to be extracted
    erpFilterType = {'bandpass'};
    erpType = {'ERP'};
    chsToPlot = {'Cz'; 'Pz'};
    
    for filt = 1 : length(erpFilterType)       
        for erp = 1 : length(erpType)      
            
            % Import the data from .MAT file
            [epochsOut, debugInfoOut] = batch_pullOut_timeDomainEpochs(fileNameFields, outlierFilenameList, erpType{erp}, erpFilterType{filt}, epochType, handles);
            
            % Preprocess the epochs for plotting
            [epochsProcessed, debugInfoOut] = batch_preProcessEpochs(epochsOut, debugInfoOut, fileNameFields, outlierFilenameList, erpType{erp}, erpFilterType{filt}, epochType, handles)
            
            % Plot
            
        end        
    end