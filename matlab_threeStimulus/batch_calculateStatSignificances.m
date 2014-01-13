function statsTests = batch_calculateStatSignificances(statsOut, normFieldName, erpBandType, erpComponent, erpFilterType, matricesSessionNorm, subjects, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempStatSignif.mat';
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

    % example input data
    noOfSessions = size(statsOut.target.dark.Fz.mean);    
    [noOfSessions, noOfSubjects] = size(statsOut.target.dark.Fz.meanIn)
    statFields = fieldnames(statsOut.target.dark.Fz)


    % construct the comparisons
    stimTypes = fieldnames(statsOut);
    conditions = fieldnames(statsOut.(stimTypes{1}));
    chNames = fieldnames(statsOut.(stimTypes{1}).(conditions{1}));

    statField = 'meanIn';

    for stim = 1 : length(stimTypes)
        
        for ch = 1 : length(chNames)
            
            for cond = 1 : length(conditions)
                comparisonMatrix(cond,:,:) = statsOut.(stimTypes{stim}).(conditions{cond}).(chNames{ch}).(statField)                
            end
            
            [statsTests.assumptions.(stimTypes{stim}).(chNames{ch}), statsTests.testResults.(stimTypes{stim}).(chNames{ch})] = ...
                stat_signifTestWrapper(comparisonMatrix, statsOut.(stimTypes{stim}), subjects, normFieldName, erpBandType, erpComponent, erpFilterType, statField, stimTypes{stim}, chNames{ch}, handles.parameters, handles)
            
        end
    end

      
    