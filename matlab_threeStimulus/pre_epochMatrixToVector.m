function EEGmat = pre_epochMatrixToVector(EEG_data, parameters, handles)

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempEpochToMat.mat';
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

   
    
    [noOfChannels, noOfSamplesPerEpoch, noOfEpochs] = size(EEG_data);
    EEGmat = zeros(noOfSamplesPerEpoch*noOfEpochs, noOfChannels);
    
    % vectorize maybe later with reshape
    for ep = 1 : noOfEpochs
        ind1 = (ep - 1)*noOfSamplesPerEpoch + 1;
        ind2 = ind1 + noOfSamplesPerEpoch - 1;
        epochPiece = (squeeze(EEG_data(:, :, ep)))';
        EEGmat(ind1:ind2, :) = epochPiece;
    end
    
    