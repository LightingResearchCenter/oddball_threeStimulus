function matrixNorm = batch_normalizeAux(matrixIn, typeOfNorm, subjects, handles)
    
    %% DEBUG
    debugMatFileName = 'tempNormAux.mat';
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
        
    if strcmp(typeOfNorm, 'darkCondition')
        normValues = batch_normalizeSubLoop(matrixIn, [], typeOfNorm, [], 'dark', 'getNorm', subjects, handles);
        matrixNorm = batch_normalizeSubLoop(matrixIn, normValues, typeOfNorm, [], 'dark', 'normalize', subjects, handles);

    elseif strcmp(typeOfNorm, 'firstSession')        
        normValues = batch_normalizeSubLoop(matrixIn, [], typeOfNorm, 1, [], 'getNorm', subjects, handles);
        matrixNorm = batch_normalizeSubLoop(matrixIn, normValues, typeOfNorm, 1, [], 'normalize', subjects, handles);
    end

    
    
    function normOut = batch_normalizeSubLoop(matrixIn, normValues, typeOfNorm, sessionToFix, conditionToFix, normPassType, subjects, handles)
        
        normOut = matrixIn;
        
        % Get field parameters        
        conditions = fieldnames(matrixIn);
        noOfConditions = length(conditions);
        auxFields = fieldnames(matrixIn.(conditions{1}));
        
        if strcmp(typeOfNorm, 'darkCondition')
            conditionFixedIndex = find(strcmp(conditions, conditionToFix) == 1);            
        elseif strcmp(typeOfNorm, 'firstSession')
            sessionFixedIndex = sessionToFix;
        else
            error('typo?')
        end
        
        for condition = 1 : noOfConditions
            
            for field = 1 : length(auxFields)
                
                auxParam = fieldnames(matrixIn.(conditions{condition}).(auxFields{field}));
                for param = 1 : length(auxParam)

                    % insideField = matrixIn.(conditions{condition}).(auxFields{field})
                    [noOfSessions, noOfSubjects] = size(matrixIn.(conditions{condition}).(auxFields{field}).(auxParam{param}));
                    
                    for session = 1 : noOfSessions

                        if strcmp(normPassType, 'getNorm')

                            if strcmp(typeOfNorm, 'darkCondition') 
                                % keep the index fixed for the
                                % condition that you save to the output
                                % cell structure
                                dataPoints = matrixIn.(conditions{conditionFixedIndex}).(auxFields{field}).(auxParam{param})(session,:);
                                normOut.(conditions{condition}).(auxFields{field}).(auxParam{param})(session,:) = dataPoints;                                                                    
                                    
                                    if session == 1 && param == 1 && field == 1 && 1 == 2
                                        disp('All the conditions should have the same values (this is the normValue)')                                        
                                        disp(['condition = "', conditions{condition}, '", session=', num2str(session), ', and the values and indices:'])
                                        disp(dataPoints)
                                    end

                            elseif strcmp(typeOfNorm, 'firstSession')        
                                
                                dataPoints = matrixIn.(conditions{condition}).(auxFields{field}).(auxParam{param})(sessionFixedIndex,:);
                                normOut.(conditions{condition}).(auxFields{field}).(auxParam{param})(session,:) = dataPoints;
                                                                
                                if condition == 1 && param == 1 && field == 1 && 1 == 2
                                    disp('All the sessions should have the same values (this is the normValue)')
                                    disp(auxParam{param})
                                    disp(auxFields{field})
                                    disp(['condition = "', conditions{condition}, '", session=', num2str(session), ', and the values and indices:'])
                                    disp(dataPoints)
                                end
                                
                            end

                        elseif strcmp(normPassType, 'normalize')

                            % does not matter know on what you have
                            % decided to normalize the data to,
                            % just matters that the first pass was
                            % correct

                            
                            
                            dataPoints_mean_in = matrixIn.(conditions{condition}).(auxFields{field}).(auxParam{param})(session,:);
                            % dataPoints_SD_in = ?;
                            norm = normValues.(conditions{condition}).(auxFields{field}).(auxParam{param})(session,:);                                            

                            [dataPoints_mean_out, ~, ~] = batch_normalizeLowLevel(dataPoints_mean_in, 'mean', [], [], norm, handles);                                         
                            % [dataPoints_SD_out, LE, UE] = batch_normalizeLowLevel(dataPoints_SD_in, 'SD', dataPoints_mean_in, dataPoints_mean_out, norm, handles);     
                            
                            % exclude outliers
                            [dataPoints_mean_out, outlierIndices] = batch_excludeOutliersDuringBatch(dataPoints_mean_out, subjects, 'normalize', handles);
                            
                            if strcmp(typeOfNorm, 'firstSession') && condition == 1 && param == 1 && field == 1 && 1 == 2
                                disp(conditions{condition})
                                disp(session)
                                disp(dataPoints_mean_in)
                                disp(norm)
                                disp(dataPoints_mean_out)
                            end

                                
                            normOut.(conditions{condition}).(auxFields{field}).(auxParam{param})(session,:) = dataPoints_mean_out;
                            %normOut.(conditions{condition}).(auxFields{field}).(auxParam{param})(session,:) = [];
                            %normOut.(conditions{condition}).(auxFields{field}).(auxParam{param})(session,:) = []; % lower error bound
                            %normOut.(conditions{condition}).(auxFields{field}).(auxParam{param})(session,:) = []; % upper error bound

                        end
                    end                    
                end
            end
        end