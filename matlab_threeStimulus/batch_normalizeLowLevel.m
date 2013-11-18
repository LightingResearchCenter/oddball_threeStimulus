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
    % subfunction, while ridiculously simple

    % to correct for "constriction effect" as the constriction causes the
    % pupil to become smaller (smaller value) then it is easier to consider
    % as amplified response when scaled to the other side of 1.    
    % onesVector = ones(1,length(in));    
    % in = onesVector - in
    % norm = onesVector - norm           
    
    % if component is negative while norm is positive, 
    % then it should be handled
    
    diff = zeros(1, length(in));
    for i = 1 : length(in)
        
        if in(i) < 0 && norm(i) > 0    
            diff(i) = in(i) - norm(i);
             % fix later
        else
            diff(i) = in(i) - norm(i);
        end
        
    end

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
        

    