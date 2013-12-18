function [heart_Out, fractalEEG_Out, eye_Out] = batch_extraSensorStats(heart, fractalEEG, eye, statsPer, stimulusType, handles)

    %% DEBUG
    debugMatFileName = 'tempStatsOutExtra.mat';
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
    
    %{
    fractalEEG.(intensity).(session).(subject).scalar
    heart.(intensity).(session).(subject).scalar
    eye.(intensity).(session).(subject).scalar
    %}
    
    % heart
    heart_Out = batch_extraProcessEachSensor(heart, 'heart', handles);
    
    % fractal EEG
    % fractalEEG_Out = batch_extraProcessEachSensor(fractalEEG, 'fractalEEG', handles)
    fractalEEG_Out = [];
    
    % eye
    % eye_Out = batch_extraProcessEachSensor(eye, 'eye', handles)
    eye_Out = [];
    
    
    function statOut = batch_extraProcessEachSensor(dataIn, sensorType, handles)
        
        conditionFields = fieldnames(dataIn);
        sessionFields = fieldnames(dataIn.(conditionFields{1}));
        subjectFields = fieldnames(dataIn.(conditionFields{1}).(sessionFields{1}));
        scalarFields = fieldnames(dataIn.(conditionFields{1}).(sessionFields{1}).(subjectFields{1}).scalar);
        vectorFields = fieldnames(dataIn.(conditionFields{1}).(sessionFields{1}).(subjectFields{1}).vector);
                
        noOfSessions = length(sessionFields);
        noOfSubjects = length(subjectFields);
        
        disp(['  .. ', sensorType])
        
        for scalarField = 1 : length(scalarFields)
        
            % Create a matrix (or a vector of the values)
            disp('     . create vector/matrix of all the subjects')
            dataOut.absolute = batch_extraCreateMatrix(dataIn, sensorType, 'scalar', handles);

            % Normalize then 
            disp('     .. normalize')
            dataOut.darkCondition = batch_extraNormalizeValues(dataOut.absolute, 'darkCondition', 'scalar', handles);
            dataOut.firstSession = batch_extraNormalizeValues(dataOut.absolute, 'firstSession', 'scalar', handles);

            % And then compute the stats of all subjects
            disp('     ... compute the stats')
            statOut = batch_extraComputeStats(dataOut, 'scalar', sensorType, handles);
            
        end