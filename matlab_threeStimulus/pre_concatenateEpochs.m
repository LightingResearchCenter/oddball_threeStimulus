function epochs_concan = pre_concatenateEpochs(epochs, parameters, handles)

    debugMatFileName = 'tempConcatenate.mat';
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
    
    % assign RTs directly to output
    epochs_concan.RT = epochs.RT;    
    
    % also the samples per epoch
    epochs_concan.samplesPerEpoch = epochs.samplesPerEpoch;
    
    
    %% PRE-ALLOCATE       

        totalLength = 0;
        [rowsIn, colsIn] = size(epochs.ERP{1});
        for i = 1 : length(epochs.ERP)
            epochs_concan.Indices{i}(1) = totalLength + 1; % start
            totalLength = totalLength + length(epochs.ERP{i});
            epochs_concan.Indices{i}(2) = totalLength; % end
        end       
        epochs_concan.ERP = zeros(totalLength,colsIn);
        
    %% ACTUALLY CONCATENATE          
        
        for i = 1 : length(epochs.ERP)
            i1 = epochs_concan.Indices{i}(1);
            i2 = epochs_concan.Indices{i}(2);
            epochs_concan.ERP(i1:i2,:) = epochs.ERP{i};
        end
        
    
        