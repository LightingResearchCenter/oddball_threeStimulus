function KSS_z = batch_zScoreTransformKSS(KSS, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempKSS_zScore.mat';
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
    
    % Group the subjects now together
    for int = 1 : length(intensityFields)    
        for sub = 1 : length(subjects)
            subjectMatrix.(subjects{sub})(:,int) = KSS.(intensityFields{int}).(subjects{sub});
        end        
    end
    
    % Z-transform each subject
    % e.g. http://www.mathworks.com/help/stats/zscore.html
    % e.g. http://www.mathworks.com/matlabcentral/newsreader/view_thread/50544
    for sub = 1 : length(subjects)
        subjectMatrix_Z.(subjects{sub}) = batch_zTransform(subjectMatrix.(subjects{sub}), handles);        
    end
        
    % Convert the subject matrix back to the input format 
    for int = 1 : length(intensityFields)    
        for sub = 1 : length(subjects)
            KSS_z.(intensityFields{int}).(subjects{sub}) = subjectMatrix_Z.(subjects{sub})(:,int);
        end        
    end
    
    
    
    %% LOW-LEVEL Z-TRANSFORM
    function matOut = batch_zTransform(matIn, handles)
        
        [sessions, intensities] = size(matIn);
        
        % vectorize
        % matIn = matIn(:)
        
        % from Statistic Toolbox
        flag = 0;
        dim = 1;
        
        try
            matOut = zscore(matIn,flag,dim);
        catch err
            % crashes if you don't have the toolbox and when someone is
            % using all the licenses
            err            
            matOut = matIn;
        end