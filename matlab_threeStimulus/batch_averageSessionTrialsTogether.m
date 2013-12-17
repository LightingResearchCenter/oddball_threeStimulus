function [matricesSessionAver, outlierOut] = batch_averageSessionTrialsTogether(matricesIntensity, fileNames, subjects, outlierFilenameList, handles)
    
    ERPtypes = fieldnames(matricesIntensity);
    conditions = fieldnames(matricesIntensity.(ERPtypes{1}));
    noOfConditions = length(conditions);
    [noOfChannels, noOfSessions, noOfTrials, noOfSubjects] = size(matricesIntensity.(ERPtypes{1}).dark);

    % matricesConcatenated = zeros(noOfChannels, noOfSessions, noOfSubjects);
    % dataMatrix = 

    dim = 1; % the return should be a row vector, one value per subject
    flag = 0;
    
    outlierOut = [];

    dataMatrixInit = 0;
    for condition = 1 : noOfConditions
        
        for j = 1 : length(ERPtypes)
            
            for ch = 1 : noOfChannels

                for session = 1 : noOfSessions                    
                                            
                    dataMatrix = zeros(length(matricesIntensity.(ERPtypes{j}).(conditions{condition})(ch, session, :, 1)), noOfSubjects);
                    dataMatrix(:,:) = NaN;

                    for subjects = 1 : noOfSubjects

                        % create the data matrix, 
                        %  -- rows : as many trials there are for session
                        %  --- cols : as many subjects that there are
                        matrixTemp = matricesIntensity.(ERPtypes{j}).(conditions{condition})(ch, session, :, subjects);
                                                
                        try
                            dataMatrix(:, subjects) = matrixTemp;
                        catch err
                            err
                            whos
                            if length(dataMatrix) < length(matrixTemp) % too many trials
                                dataMatrix(:, subjects) = matrixTemp(1,1,1:length(dataMatrix));
                                warning('Too many trials per ERP TYPE, bug at the pre_findERP_Epochs() or in the original experiment. I.e. you should always have 40 trials per target, and the artifacts are just NaNs, but the length should be the same')
                            elseif length(dataMatrix) < length(matrixTemp) % too little trials
                                dataMatrix(1:length(matrixTemp), subjects) = matrixTemp(1,1,:);
                                warning('Too little trials per ERP TYPE, bug at the pre_findERP_Epochs() or in the original experiment. I.e. you should always have 40 trials per target, and the artifacts are just NaNs, but the length should be the same')
                            else
                                error('Why are we here actually? Sizes are the same yes')
                            end
                            
                        end

                    end 

                    % calculate the stats now from this 2D Matrix, for each
                    % channel independently
                    % dataMatrix;
                    % av1 = batch_calculateStatsForMatrix(dataMatrix, dim, flag, handles);
                    % ERPtypes{j}
                    dataMatrix = batch_excludeOutliersBeforeStats(dataMatrix, subjects, handles)
                    matricesSessionAver.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){session} = batch_calculateStatsForMatrix(dataMatrix, dim, flag, handles);

                    % session

                end              


            end
        end
    end
