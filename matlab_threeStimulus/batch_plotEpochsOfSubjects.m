function batch_plotEpochsOfSubjects(epochsMatrix, erpType, erpFilterType, epochType, fileNameFields, statParam, handles)

    %% DEBUG
    debugMatFileName = 'tempBatchEpochPlots.mat';
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
    parameters = handles.parameters;

    
    % get fieldnames IN
    conditions = fieldnames(epochsMatrix)
    sessions = fieldnames(epochsMatrix.(conditions{1}))
    erpResponses = fieldnames(epochsMatrix.(conditions{1}).(sessions{1}))
    erpTypes = fieldnames(epochsMatrix.(conditions{1}).(sessions{1}).(erpResponses{1}))
    filterTypes = fieldnames(epochsMatrix.(conditions{1}).(sessions{1}).(erpResponses{1}).(erpTypes{1}))
    statFields = fieldnames(epochsMatrix.(conditions{1}).(sessions{1}).(erpResponses{1}).(erpTypes{1}).(filterTypes{1}))
    
    exIn = epochsMatrix.(conditions{1}).(sessions{1}).(erpResponses{1}).(erpTypes{1}).(filterTypes{1})
    
    [noOfDataSamplesPerERP, noOfChannels, noOfSubjects] = size(epochsMatrix.(conditions{1}).(sessions{1}).(erpResponses{1}).(erpTypes{1}).(filterTypes{1}).mean)
    t = 1000*(linspace(-parameters.oddballTask.ERP_baseline, parameters.oddballTask.ERP_duration, noOfDataSamplesPerERP))';
        
    
    
    %% SUBJECTS COMPARED with AVERAGES
    %{
    batch_plot_subjectERPs_withAverage(conditions, sessions, erpResponses, erpTypes, filterTypes, statFields, statParam, ...
                                       noOfDataSamplesPerERP, noOfChannels, noOfSubjects, ...
                                       t, epochsMatrix, parameters, handles)    
    %}
        
        
    %% CHANNELS COMPARED (AVERAGE)
    rowParameter = conditions;
    colParameter = 1 : parameters.EEG.nrOfChannels;
    gcaParameter = erpResponses;
    
    batch_plot_ERPs_perSomeParameter(rowParameter, colParameter, gcaParameter, ...
                                     conditions, sessions, erpResponses, erpTypes, filterTypes, statFields, statParam, ...
                                       noOfDataSamplesPerERP, noOfChannels, noOfSubjects, ...
                                       t, epochsMatrix, parameters, handles)
    
    
    
        
    %% COMPARE THE CONDITIONS (fix this?)
    rowParameter = erpResponses;
    colParameter = 1 : parameters.EEG.nrOfChannels;
    gcaParameter = conditions;
    
    batch_plot_ERPs_perSomeParameter(rowParameter, colParameter, gcaParameter, ...
                                     conditions, sessions, erpResponses, erpTypes, filterTypes, statFields, statParam, ...
                                       noOfDataSamplesPerERP, noOfChannels, noOfSubjects, ...
                                       t, epochsMatrix, parameters, handles)
    