function batchCompare_KSS(handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempKSS.mat';
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

    % IMPORT
    KSS = batch_importKSS(handles.path.dataKSS, handles);
    
    % Normalize (z-score), within subjects
    KSS_z = batch_zScoreTransformKSS(KSS, handles);
    
    % Normalize to the first session    
    KSS_Norm = batch_KSS_normalize(KSS, handles);
    KSS_z_Norm = batch_KSS_normalize(KSS_z, handles);    
    
    % PLOT
    batch_plotKSS(KSS, KSS_z, KSS_Norm, KSS_z_Norm, handles.parameters, handles)