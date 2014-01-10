function [stat_assumptions, stat_results] = stat_signifTestWrapper(inputMatrix, subjects, normFieldName, stimType, chName, parameters, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempStatSignifWrapper.mat';
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

    stimType
    chName

    [noOfConditions, noOfSessions, noOfSubjects] = size(inputMatrix);
    
    noOfSubgroups = noOfConditions * noOfSessions;
    
    parameters.stats.shapWilk_pThreshold = 0.05;
    parameters.stats.bartlett_pThreshold = 0.05;
    
    conditions = {'dark'; 'dim'; 'bright'};
    sessions = {'session1'; 'session2'; 'session3'; 'session4'};
        
   
    %% TEST ASSUMPTIONS
    % 1) Normal distribution
    % 2) Homogeneity of variances
    % 3) Sphericity 
    for condition = 1 : noOfConditions       
        for session = 1 : noOfSessions           
            stat_assumptions{condition, session} = batch_statTestAssumptions(squeeze(inputMatrix(condition, session, :)), inputMatrix, conditions, condition, session, parameters, handles);
        end
    end
    

    %% Actual TESTS   
    stat_results = batch_statTestResults(inputMatrix, stat_assumptions, normFieldName, stimType, chName, conditions, sessions, subjects, parameters, handles);

    