function [matrixOut, matrixRegress, filtOnly, NaN_indices] = pre_componentArtifactFiltering(matrixIn, rawMax, rowsIn, colsIn, loFreq, hiFreq, filterOrder, filterOrderSteep, dataType, parameters, handles)   

    
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempComponentArtifact.mat';
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
    

    %% Notch filter the mains (60 Hz) 
    % out using the modified filter from BioSig
    
        if strcmp(dataType, 'General') == 1
            
            matrixOut1a = zeros(rowsIn, colsIn);  % e.g. a lot of rows x 6 channels            
            HDR.SampleRate = parameters.EEG.srate;

                disp(['   remove 50/60 Hz (Notch, from BioSig package)'])
                for j = 1 : colsIn                        
                    matrixOut1a(:,j) = pre_remove5060hz_modified(matrixIn(:,j), HDR, 'PCA 60');                
                end
                dataIn = matrixOut1a; % for debug plotting
                filtOnly = matrixOut1a; 

        else
            matrixOut1a = matrixIn;                    
        end
        
    %% BANDPASS FILTER
    
        [rowsIn, colsIn] = size(matrixOut1a);
        matrixOut1 = zeros(rowsIn, colsIn);  % e.g. a lot of rows x 6 channels            
        
            disp(['    - BANDPASS FILTERING'])
            disp(['        .. ', dataType,  ' (', num2str(loFreq), '-', num2str(hiFreq), ' Hz), order = ', num2str(filterOrder), ', ZeroPhase'])
            for j = 1 : colsIn                       
                matrixOut1(:,j) = pre_bandbassFilter(matrixOut1a(:,j), parameters.EEG.srate, [hiFreq, loFreq], filterOrder, filterOrderSteep, handles);                
            end        
    
            filtOnly = matrixOut1; 
    
        
    %% FASTER for ARTIFACT REJECTION & CORRECTION
    
        % Use easier variable names        
        EEG = matrixOut1(:,1:handles.parameters.EEG.nrOfChannels);
        EOG = matrixOut1(:,handles.parameters.EEG.nrOfChannels+1);
        ECG = matrixOut1(:,handles.parameters.EEG.nrOfChannels+2);
        
    if parameters.artifacts.continuousFASTER == 1 && strcmp(dataType, 'General')
        
        parameters.artifacts.contFaster_regressEOGbeforeFaster = 1;
        
        %% NOTE! the epoching trims a bit the data to integer number of epochLenghts
        
        % Rethink a bit the procedure, and make it a better for continous
        % data after Christmas, does not leave that much for the
        % epoch-by-epoch Faster to work with
        
        
        % If you wanna do regress EOG before whole FASTER
        if parameters.artifacts.contFaster_regressEOGbeforeFaster == 1
            hiFreqEOG = parameters.artifacts.regressBirectContamCutoffs(2);
            loFreqEOG = parameters.artifacts.regressBirectContamCutoffs(1);      
            filterOrder = parameters.filterOrder;        
            disp(['            - [beforeFASTER] mitigate bidirectional contamination with regressEOG by filtering EOG (', num2str(loFreqEOG), '-', num2str(hiFreqEOG), ' Hz) before regression'])            
            EOG = pre_bandbassFilter(EOG, parameters.EEG.srate, [hiFreqEOG, loFreqEOG], filterOrder, [], handles); 
            disp(['             ... EOG REGRESSION'])
            EEG = pre_artifactRegressionWrapper(EEG, EOG, parameters.EEG.nrOfChannels, parameters);
        end
        
        % Now split the continuous data into chunks and evaluate each chunk
        % as an epoch, and determine whether this chunk is artifacted or
        % not
        [matrixOut, ~, NaN_indices] = pre_FASTER_forContinuousData(EEG, EOG, ECG, parameters, handles);
        matrixRegress = [EEG EOG ECG];
        
        disp(['          ARTIFACT SUMMARY:'])
        for ch = 1: parameters.EEG.nrOfChannels
            
            isNanSum = sum(NaN_indices(:,ch));
            isNanRatio = isNanSum / length(NaN_indices(:,ch));
            disp(['           - ', parameters.BioSemi.chName{ch+parameters.BioSemi.chOffset}, ': ' ...
                    num2str(100*isNanRatio,4), '% of the input is artifacts'])
        end        
        
    elseif strcmp(dataType, 'General')
        
        hiFreqEOG = parameters.artifacts.regressBirectContamCutoffs(2);
        loFreqEOG = parameters.artifacts.regressBirectContamCutoffs(1);      
        filterOrder = parameters.filterOrder;        
        disp(['            - [beforeFASTER] mitigate bidirectional contamination with regressEOG by filtering EOG (', num2str(loFreqEOG), '-', num2str(hiFreqEOG), ' Hz) before regression'])            
        EOG = pre_bandbassFilter(EOG, parameters.EEG.srate, [hiFreqEOG, loFreqEOG], filterOrder, [], handles); 
        disp(['             ... EOG REGRESSION'])
        EEG = pre_artifactRegressionWrapper(EEG, EOG, parameters.EEG.nrOfChannels, parameters);
        
        matrixRegress = [EEG EOG ECG];
        matrixOut = matrixRegress;
        NaN_indices = false(length(matrixOut), parameters.EEG.nrOfChannels);
       
    else
        
        matrixOut = matrixOut1;
        matrixRegress = matrixOut;
        NaN_indices = false(length(matrixOut), parameters.EEG.nrOfChannels);
        
    end
        
        
    