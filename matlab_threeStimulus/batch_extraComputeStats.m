function statsOut = batch_extraComputeStats(dataIn, scalarOrVector, sensorType, handles)

    %% DEBUG
    debugMatFileName = 'tempStatExtra.mat';
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

    normType = fieldnames(dataIn);
    conditionFields = fieldnames(dataIn.(normType{1}));
    dataTypes = fieldnames(dataIn.(normType{1}).(conditionFields{1}).(scalarOrVector));
    
    % Is cell for EEG data where each cell is one EEG channel, whereas
    % there are no cells for heart and eye movement data
    if iscell(dataIn.(normType{1}).(conditionFields{1}).(scalarOrVector).(dataTypes{1}))
        dataSamples = dataIn.(normType{1}).(conditionFields{1}).(scalarOrVector).(dataTypes{1}){1};
        [noOfDataPointsPerPersonPerSession, noOfSessions, noOfSubjects] = size(dataSamples);
        isCell = 1;
    else
        dataSamples = dataIn.(normType{1}).(conditionFields{1}).(scalarOrVector).(dataTypes{1});
        [noOfDataPointsPerPersonPerSession, noOfSessions, noOfSubjects] = size(dataSamples);
        isCell = 0;
    end
        
    % Now the input data is 3D and we need to "get rid" of the last
    % dimension which is the subjects so we would average over the subjects    
    for norm = 1 : length(normType)
       
        for condition = 1 : length(conditionFields)
            
            for type = 1 : length(dataTypes)
                
                % for session = 1 : noOfSessions
                
                    % for EEG (channels)
                    if isCell == 1                       
                        
                        for ch = 1 : length(dataIn.(normType{norm}).(conditionFields{condition}).(scalarOrVector).(dataTypes{type}))
                            dataSamplesPerSession = dataIn.(normType{norm}).(conditionFields{condition}).(scalarOrVector).(dataTypes{type}){ch}(:,:,:);
                            dataSamplesPerSession = squeeze(dataSamplesPerSession); % squeeze and remove the singleton (session)
                            if noOfDataPointsPerPersonPerSession > 1
                                statsOut.(normType{norm}).(conditionFields{condition}).(scalarOrVector).(dataTypes{type}) = batch_calculateStatsForMatrix(dataSamplesPerSession, 2, 0, handles);
                            else
                                statsOut.(normType{norm}).(conditionFields{condition}).(scalarOrVector).(dataTypes{type}) = batch_calculateStatsForVector(dataSamplesPerSession, 1, handles);
                            end
                        end
                        
                    % For others
                    else
                        
                        dataSamplesPerSession = dataIn.(normType{norm}).(conditionFields{condition}).(scalarOrVector).(dataTypes{type})(:,:,:);
                        dataSamplesPerSession = squeeze(dataSamplesPerSession); % squeeze and remove the singleton (session)
                        
                        if noOfDataPointsPerPersonPerSession > 1
                            statsOut.(normType{norm}).(conditionFields{condition}).(scalarOrVector).(dataTypes{type}) = batch_calculateStatsForMatrix(dataSamplesPerSession, 2, 0, handles);
                        else
                            statsOut.(normType{norm}).(conditionFields{condition}).(scalarOrVector).(dataTypes{type}) = batch_calculateStatsForVector(dataSamplesPerSession, 1, handles);
                        end
                        
                    end    
                    
                    % pause
                    
                % end                
                
            end
            
        end
        
    end
    