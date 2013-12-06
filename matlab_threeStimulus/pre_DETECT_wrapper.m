function indicesArtifact = pre_DETECT_wrapper(EEG, EOG, ECG, j, dataType, parameters, handles)
    
    debugMatFileName = 'tempDETECTwrapper.mat';
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
    fdy = 1