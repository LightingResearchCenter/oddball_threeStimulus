function epochs = pre_deconcatenateEpochs(epochs_concan, parameters, handles)

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if 2 == 1 % handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempDeconcenateEpochs.mat';
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
    
    % assign RTs directly to output
    epochs.RT = epochs_concan.RT;
    
    % also the samples per epoch
    epochs.samplesPerEpoch = epochs_concan.samplesPerEpoch;
    
    % same for indices
    epochs.Indices = epochs_concan.Indices;
            
    % DeCONCATENATE
        
        for i = 1 : length(epochs_concan.Indices)
            i1 = epochs_concan.Indices{i}(1);
            i2 = epochs_concan.Indices{i}(2);
            epochs.ERP{i} = epochs_concan.ERP(i1:i2,:);
        end        
        