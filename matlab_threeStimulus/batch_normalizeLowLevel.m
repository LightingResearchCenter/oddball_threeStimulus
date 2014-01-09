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
    
    % if you are dealing with negative values, the above formulation will
    % give you always suppression, maybe make this a bit more elegant later
    % e.g. when input is 0.5 and norm is -1.0, then diff = 0.5 - (-1.0) =
    % 1.5, and out is 1.5 / -1.0 = -1.5, when the change then is actually
    % positive
    
    for i = 1 : length(norm)
       
        negBoolean(i,1) = logical(in(i) > norm(i) && norm(i) < 0);
        if negBoolean(i,1) == 1
            negMultip(i,1) = -1;
        else
            negMultip(i,1) = 1;
        end         
        out(i) = negMultip(i) * out(i);
        
    end
    
    %[in diff norm out negBoolean negMultip]
    %pause
    
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
        

    