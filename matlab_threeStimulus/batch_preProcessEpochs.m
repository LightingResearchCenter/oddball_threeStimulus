function [epochsMatrix, debugInfoOut] = batch_preProcessEpochs(epochsStat, debugInfoOut, fileNameFields, outlierFilenameList, erpType, erpFilterType, epochType, handles)

    %% DEBUG
    debugMatFileName = 'tempPreprocessTimeDomainEpochs.mat';
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
    
    whos
    exInput = epochsStat.dim.session1.target.(erpType).ka.bandpass
    
    %{
    % exInput = 

      mean: [5734x4 double]
    median: [5734x4 double]
        SD: [5734x4 double]
         n: [30 31 23 31]
    %}
            
    % create matrices from structure monsters (i.e. get rid of the subject
    % field name), e.g.    
    epochsMatrix = batch_epochStructToMatrix(epochsStat);
    
        % plot subjects
        statParam = 'mean';
        batch_plotEpochsOfSubjects(epochsMatrix, erpType, erpFilterType, epochType, fileNameFields, statParam, handles)
    
    % average the subjects
    %epochsAveraged = 
    
    
    
    
    function epochsMatrix = batch_epochStructToMatrix(epochsStat)
        
        % get fieldnames IN
        conditions = fieldnames(epochsStat);
        sessions = fieldnames(epochsStat.(conditions{1}));
        erpResponses = fieldnames(epochsStat.(conditions{1}).(sessions{1}));
        erpTypes = fieldnames(epochsStat.(conditions{1}).(sessions{1}).(erpResponses{1}));
        subjects = fieldnames(epochsStat.(conditions{1}).(sessions{1}).(erpResponses{1}).(erpTypes{1}));
        filterTypes = fieldnames(epochsStat.(conditions{1}).(sessions{1}).(erpResponses{1}).(erpTypes{1}).(subjects{1}));
        statFields = fieldnames(epochsStat.(conditions{1}).(sessions{1}).(erpResponses{1}).(erpTypes{1}).(subjects{1}).(filterTypes{1}));
        
        for cond = 1 : length(conditions)
            for ses = 1 : length(sessions)
                for erpResp = 1 : length(erpResponses)
                    for erp = 1 : length(erpTypes)
                        for filt = 1 : length(filterTypes)
                            for stat = 1 : length(statFields)
                                for sub = 1 : length(subjects)                                       
                                    statMatrix = epochsStat.(conditions{cond}).(sessions{ses}).(erpResponses{erpResp}).(erpTypes{erp}).(subjects{sub}).(filterTypes{filt}).(statFields{stat});
                                    epochsMatrix.(conditions{cond}).(sessions{ses}).(erpResponses{erpResp}).(erpTypes{erp}).(filterTypes{filt}).(statFields{stat})(:,:,sub) = statMatrix;
                                end                                
                            end
                        end
                    end
                end
            end
        end
            