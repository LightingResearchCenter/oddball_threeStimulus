function [epochsProcessed, debugInfoOut] = batch_preProcessEpochs(epochsStat, debugInfoOut, fileNameFields, outlierFilenameList, erpType, erpFilterType, epochType, handles)

    %% DEBUG
    debugMatFileName = 'tempPreprocessTimeDomainEpochs.mat';
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
    
    epochsProcessed = []
       