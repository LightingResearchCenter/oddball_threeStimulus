function [out, LE, UE] = batch_normalizeLowLevel(in, dataType, meanIn, meanOut, norm, handles)
    
    %% DEBUG
    %{
    debugMatFileName = 'tempNormLowLevel.mat';
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
    %}

    % if you wanna play with this as well, so that is why this is a
    % subfunction, while rather simple
    
    %diff = zeros(1, length(in));
    diff = in - norm;
    out = diff ./ norm;
    
    if strcmp(dataType, 'SD')
        %dataType
        %meanIn
        %meanOut
        LE = [];
        UE = [];
    else
        LE = [];
        UE = [];
    end
        

    