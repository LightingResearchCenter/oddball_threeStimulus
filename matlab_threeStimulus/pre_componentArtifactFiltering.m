function [matrixOut, filtOnly, NaN_indices] = pre_componentArtifactFiltering(matrixIn, artifactIndices, rawMax, rowsIn, colsIn, loFreq, hiFreq, filterOrder, filterOrderSteep, dataType, parameters, handles)   

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

    %% Notch filter the mains (60 Hz) 
    % out using the modified fitler from BioSig
    if strcmp(dataType, 'General') == 1
        matrixOut1a = zeros(rowsIn, colsIn);  % e.g. a lot of rows x 6 channels            
        HDR.SampleRate = parameters.EEG.srate;

            disp(['   remove 50/60 Hz (Notch, from BioSig package)'])
            for j = 1 : colsIn                        
                matrixOut1a(:,j) = pre_remove5060hz_modified(matrixIn(:,j), HDR, 'PCA 60');                
            end
    else
        matrixOut1a = matrixIn;
    end
    
        
    %% Remove here fixed threshold artifacts (or get only the indices, 
    % and not pass the NaN values to the bandpass filter
    if strcmp(dataType, 'General') == 1
        [~, NaN_indices, numberOfNaNs] = pre_artifactFixedThreshold(matrixOut1a(:,1:handles.parameters.EEG.nrOfChannels), matrixOut1a(:,parameters.EEG.nrOfChannels + 1), parameters, handles);
    else
        NaN_indices = artifactIndices;
    end
    
    %% band-pass filter the data using a a subfunction            
    matrixOut1 = zeros(rowsIn, colsIn);  % e.g. a lot of rows x 6 channels            
    
        for j = 1 : colsIn                        
            matrixOut1(:,j) = pre_bandbassFilter(matrixOut1a(:,j), parameters.EEG.srate, [hiFreq, loFreq], filterOrder, filterOrderSteep, handles);                
        end
        
        filtOnly = matrixOut1;
        
    %% remove artifacts (i.e. deviant voltages, ECG and EOG artifacts)
    if parameters.artifacts.applyRegressEOG == 1
        
        disp(['   Remove the artifacts (', dataType, ') threshold: ', num2str(parameters.artifacts.fixedThr), ' uV, bandpass: ', num2str(loFreq), '-', num2str(hiFreq), ' Hz'])                
        
        % correct the NanIndices to have 6 columns
        [rows,cols] = size(NaN_indices);
        NaN_indices_corr = zeros(rows,cols+2);
        NaN_indices_corr(:, 1:parameters.EEG.nrOfChannels) = NaN_indices;
        NaN_indices_corr = logical(NaN_indices_corr);
        
        % remove now the fixed threshold artifacts and do the other
        % operations
        matrixOut1(NaN_indices_corr) = NaN;        
        
        % Use easier variable names
        EEG = matrixOut1(:,1:parameters.EEG.nrOfChannels);
        EOG = matrixOut1(:,parameters.EEG.nrOfChannels + 1);
        ECG = matrixOut1(:,parameters.EEG.nrOfChannels + 2);
        
        % use subfunction
        matrixOut = matrixOut1;
        matrixOut(:,1:parameters.EEG.nrOfChannels) = pre_artifactRemovalInputData(EEG, EOG, ECG, j, dataType, parameters, handles);
                                                                       
        % debug info
        filtMax = max(max(abs(matrixOut(:,1:parameters.EEG.nrOfChannels))));
        disp(['       .. Max voltage of "', dataType, '" matrix: ', num2str(filtMax), ' uV, <-- from ', num2str(rawMax), ' uV'])
                
    else
        
        disp(['   ', dataType, ': Not applying the regress_eog() -based artifact removal'])
        
    end
    
   
    