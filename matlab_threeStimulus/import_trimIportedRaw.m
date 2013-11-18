function [dataMatrix, triggers] = import_trimIportedRaw(dataMatrix, triggers)

    debugMatFileName = 'tempTrim.mat';
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
    %triggers(1:100)

    unique(triggers)
    % plot(triggers)