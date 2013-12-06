function [matrixOut, filtOnly, NaN_indices] = pre_componentArtifactFiltering(matrixIn, artifactIndices, rawMax, rowsIn, colsIn, loFreq, hiFreq, filterOrder, filterOrderSteep, dataType, parameters, handles)   

    %{
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
    %}

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
        % noOfNans1 = sum(sum(isnan(matrixOut1a(:,1:handles.parameters.EEG.nrOfChannels)))) % there shouldn't be any at this point
    else
        NaN_indices = artifactIndices;
    end
    
    %% band-pass filter the data using a a subfunction            
    matrixOut1 = zeros(rowsIn, colsIn);  % e.g. a lot of rows x 6 channels            
    
        for j = 1 : colsIn                        
            matrixOut1(:,j) = pre_bandbassFilter(matrixOut1a(:,j), parameters.EEG.srate, [hiFreq, loFreq], filterOrder, filterOrderSteep, handles);                
        end
        
        filtOnly = matrixOut1;
        
    % if parameters.artifacts.useFASTER == 1
        % we remove the artifacts for epochs later using FASTER, see the call from
        % PROCESS_singleFile()
        matrixOut = matrixOut1;

        % disp(['   Skipping BioSig Artifact removal, using FASTER later (', dataType, ') threshold: ', num2str(parameters.artifacts.fixedThr), ' uV, bandpass: ', num2str(loFreq), '-', num2str(hiFreq), ' Hz'])                

        % NOTE, now we have two artifact detection pipelines at the same
        % really, as the following BioSig-based artifact correction /
        % removal are used for time series EEG analysis (power analysis,
        % multifractality analysis, etc. what you can think of for the full
        % EEG recording

        % And then we feed the non-artifact corrected (only bandpass + notch
        % filtered) epochs to FASTER artifact correction which then outputs
        % those epochs for time-frequency analysis, and for the ERP
        % component analysis

        %% YOU COULD OPTIMIZE LATER THIS!!! REMOVE REDUNDANCY        
        
   
    %% remove artifacts (i.e. deviant voltages, ECG and EOG artifacts)
    
        % correct the NanIndices to have 6 columns
            [rows,cols] = size(NaN_indices);
            NaN_indices_corr = zeros(rows,cols+2);
            NaN_indices_corr(:, 1:parameters.EEG.nrOfChannels) = NaN_indices;
            NaN_indices_corr = logical(NaN_indices_corr);
    
        % Use easier variable names
            EEG = matrixOut1(:,1:parameters.EEG.nrOfChannels);
            EOG = matrixOut1(:,parameters.EEG.nrOfChannels + 1);
            ECG = matrixOut1(:,parameters.EEG.nrOfChannels + 2);
    
        if parameters.artifacts.useDETECT == 1
            
            disp(['   Remove the artifacts with DETECT'])            
            indicesArtifact = pre_DETECT_wrapper(EEG, EOG, ECG, j, dataType, parameters, handles);
    
        elseif parameters.artifacts.applyRegressEOG == 1

            disp(['   Remove the artifacts (', dataType, ') threshold: ', num2str(parameters.artifacts.fixedThr), ' uV, bandpass: ', num2str(loFreq), '-', num2str(hiFreq), ' Hz'])                                 

            % use subfunction
            matrixOut = matrixOut1;
            matrixOut(:,1:parameters.EEG.nrOfChannels) = pre_artifactRemovalInputData(EEG, EOG, ECG, j, dataType, parameters, handles);

            % debug info
            filtMax = max(max(abs(matrixOut(:,1:parameters.EEG.nrOfChannels))));
            disp(['       .. Max voltage of "', dataType, '" matrix: ', num2str(filtMax), ' uV, <-- from ', num2str(rawMax), ' uV'])

        else

            disp(['   ', dataType, ': Not applying the regress_eog() -based artifact removal'])

        end

    if strcmp(dataType, 'General') == 1
        % Don't convert the fixed threshold artifacts to NaNs
        matrixOut = matrixOut1;

    else
        % remove now the fixed threshold artifacts and do the other
        % operations
        matrixOut(NaN_indices_corr) = NaN;
    end
    