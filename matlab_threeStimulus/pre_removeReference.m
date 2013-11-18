function dataMatrix = pre_removeReference(refRightEar, refLeftEar, dataMatrixIn, offset, parameters, handles)

    debugMatFileName = 'tempReferenceRemoval.mat';
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

    % add some check here at some point for the data quality, as the other
    % channel could be quite bad if the electrode starts to slip

    ref = (refRightEar + refLeftEar) / 2;
    [rowsIn, colsIn] = size(dataMatrixIn); % e.g. a lot of rows x 8 channels
    dataMatrix = zeros(rowsIn, colsIn-offset);

    for j = 1+offset : colsIn
        dataMatrix(:,j-offset) = dataMatrixIn(:,j) - ref;
    end
    
    % DEBUG
    if handles.flags.showDebugPlots == 1
        plot_debugReferenceChannels(ref, refRightEar, refLeftEar, parameters, handles)
    end
    