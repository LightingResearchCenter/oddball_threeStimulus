function KSS_norm = batch_KSS_normalize(KSS, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempKSSnormalize.mat';
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
    
    % get data fields
    intensityFields = fieldnames(KSS);
    subjects = fieldnames(KSS.(intensityFields{1}));       
    
    % normalize to the first session
    sessionIndex = 1;
    for int = 1 : length(intensityFields)    
        for sub = 1 : length(subjects)            
            noOfSessions = length(KSS.(intensityFields{int}).(subjects{sub}));
            for session = 1 : noOfSessions
                normMatrix.(intensityFields{int}).(subjects{sub})(session,1) = KSS.(intensityFields{int}).(subjects{sub})(sessionIndex);
            end
            
            in = KSS.(intensityFields{int}).(subjects{sub});
            norm = normMatrix.(intensityFields{int}).(subjects{sub});
            [KSS_norm.(intensityFields{int}).(subjects{sub}), ~, ~] = batch_normalizeLowLevel(in, 'KSS', [], [], norm, handles);
            
            %out = KSS_norm.(intensityFields{int}).(subjects{sub});            
            %pause
            
        end        
    end