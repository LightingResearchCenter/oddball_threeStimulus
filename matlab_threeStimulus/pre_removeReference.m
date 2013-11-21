function dataMatrix = pre_removeReference(refRightEar, refLeftEar, dataMatrixIn, offset, parameters, handles)

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
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
    end

    % add some check here at some point for the data quality, as the other
    % channel could be quite bad if the electrode starts to slip

    % additionally there could be spurious artifact only for one ear if the
    % subject for example corrects the other headphone during the session

    % Remove the power line interference    
    HDR.SampleRate = parameters.EEG.srate;    
    ref_right_notch60 = pre_remove5060hz_modified(refRightEar, HDR, 'NOTCH 60');
    ref_left_notch60 = pre_remove5060hz_modified(refLeftEar, HDR, 'NOTCH 60');
    
        if handles.flags.showDebugPlots == 1
            figure('Color', 'w')
            hold on
            plot(refRightEar, 'r')
            plot(ref_right_notch60, 'b')
            hold off
            legend('In', 'After Notch60')
        end

    % add denoising filter if you want to have more of a general trend?
    % Wavelet denoising or something, bandpass filter would remove the DC
    % offset that we want to keep
    ref = (ref_right_notch60 + ref_left_notch60) / 2;
    [rowsIn, colsIn] = size(dataMatrixIn); % e.g. a lot of rows x 8 channels
    dataMatrix = zeros(rowsIn, colsIn-offset);

    for j = 1+offset : colsIn
        dataMatrix(:,j-offset) = dataMatrixIn(:,j) - ref;
    end
    
    % DEBUG
    if handles.flags.showDebugPlots == 1
        plot_debugReferenceChannels(ref, refRightEar, refLeftEar, parameters, handles)
    end
    