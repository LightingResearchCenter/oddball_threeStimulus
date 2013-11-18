function matrixNorm = batch_extraNormalizeValues(matrixIn, typeOfNorm, scalarOrVector, handles)
    
    %% DEBUG
    debugMatFileName = 'tempNormExtra.mat';
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
        normValues = batch_extraNormalizeSubLoop(matrixIn, [], typeOfNorm, [], 'dark', 'getNorm', handles);
        matrixNorm = batch_extraNormalizeSubLoop(matrixIn, normValues, typeOfNorm, [], 'dark', 'normalize', handles);

    elseif strcmp(typeOfNorm, 'firstSession')        
        normValues = batch_extraNormalizeSubLoop(matrixIn, [], typeOfNorm, 1, [], 'getNorm', handles);
        matrixNorm = batch_extraNormalizeSubLoop(matrixIn, normValues, typeOfNorm, 1, [], 'normalize', handles);
    end

    
    
    function normOut = batch_extraNormalizeSubLoop(matrixIn, normValues, typeOfNorm, sessionToFix, conditionToFix, normPassType, handles)
        
        normOut = matrixIn;
        
        % Get field parameters        
        conditionFields = fieldnames(matrixIn);
        noOfConditions = length(conditionFields);
        dataTypes = fieldnames(matrixIn.(conditionFields{1}));
        varFields = fieldnames(matrixIn.(conditionFields{1}).(dataTypes{1}));
        
        if strcmp(typeOfNorm, 'darkCondition')
            conditionFixedIndex = find(strcmp(conditionFields, conditionToFix) == 1);            
        elseif strcmp(typeOfNorm, 'firstSession')
            sessionFixedIndex = sessionToFix;
        else
            error('typo?')
        end
        
        for condition = 1 : noOfConditions % dark / dim / bright
            
            for field = 1 : length(dataTypes) % scalar / vector                               
                    
                variableFields = fieldnames(matrixIn.(conditionFields{condition}).(dataTypes{field}));
                
                for variable =  1 : length(variableFields)
                                
                    isCellOrNot = iscell(matrixIn.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable}));
                    
                    if isCellOrNot == 1
                        %isCellOrNot
                        [noOfDataPoints, noOfSessions, noOfSubjects] = size(matrixIn.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable}){1});
                    else
                        %isCellOrNot
                        [noOfDataPoints, noOfSessions, noOfSubjects] = size(matrixIn.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable}));
                    end

                    for session = 1 : noOfSessions

                        if strcmp(normPassType, 'getNorm')

                            if strcmp(typeOfNorm, 'darkCondition') 
                                % keep the index fixed for the
                                % condition that you save to the output
                                % cell structure
                                if isCellOrNot == 1
                                    for ch = 1 : length(matrixIn.(conditionFields{conditionFixedIndex}).(dataTypes{field}).(variableFields{variable}))
                                        dataPoints = matrixIn.(conditionFields{conditionFixedIndex}).(dataTypes{field}).(variableFields{variable}){ch}(:,session,:);
                                        normOut.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable}){ch}(:,session,:) = dataPoints;
                                    end
                                else
                                    dataPoints = matrixIn.(conditionFields{conditionFixedIndex}).(dataTypes{field}).(variableFields{variable})(:,session,:);
                                    normOut.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable})(:,session,:) = dataPoints;    
                                end

                                    if session == noOfSessions && variable == 1 && field == 1
                                        disp('All the conditionFields should have the same values (this is the normValue)')                                        
                                        disp(['condition = "', conditionFields{condition}, '", session=', num2str(session), ', and the values and indices:'])
                                        disp(normOut.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable})(:,session,:))                                        
                                        disp(squeeze(dataPoints))
                                    end

                            elseif strcmp(typeOfNorm, 'firstSession')        

                                if isCellOrNot == 1
                                    for ch = 1 : length(matrixIn.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable}))
                                        dataPoints = matrixIn.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable}){ch}(:,sessionFixedIndex,:);    
                                        normOut.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable}){ch}(:,session,:) = dataPoints;
                                    end
                                else
                                    dataPoints = matrixIn.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable})(:,sessionFixedIndex,:);
                                    normOut.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable})(:,session,:) = dataPoints;
                                end

                                %dataPoints = matrixIn.(conditions{condition}).(auxFields{field}).(auxParam{param})(sessionFixedIndex,:);
                                %normOut.(conditions{condition}).(auxFields{field}).(auxParam{param})(session,:) = dataPoints;
                                
                                if condition == 1 && variable == 1 && field == 1 
                                    disp('All the sessions should have the same values (this is the normValue)')
                                    disp(variableFields{variable})
                                    disp(dataTypes{field})
                                    disp(['condition = "', conditionFields{condition}, '", session=', num2str(session), ', and the values and indices:'])
                                    disp(squeeze(dataPoints))
                                end

                            end

                        elseif strcmp(normPassType, 'normalize')

                            % does not matter know on what you have
                            % decided to normalize the data to,
                            % just matters that the first pass was
                            % correct
                            
                            if isCellOrNot == 1
                                
                                for ch = 1 : length(matrixIn.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable}))
                                    dataPoints_mean_in = matrixIn.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable}){ch}(:,session,:);                            
                                    norm = normValues.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable}){ch}(:,session,:);
                                    dataPoints_mean_out = dataPoints_mean_in ./ norm;
                                    normOut.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable}){ch}(:,session,:) = dataPoints_mean_out;
                                end
                                
                            else
                                dataPoints_mean_in = matrixIn.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable})(:,session,:);                            
                                norm = normValues.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable})(:,session,:);
                                dataPoints_mean_out = dataPoints_mean_in ./ norm;
                                normOut.(conditionFields{condition}).(dataTypes{field}).(variableFields{variable})(:,session,:) = dataPoints_mean_out;
                            end

                            if strcmp(typeOfNorm, 'firstSession') && condition == 1 && variable == 1 && field == 1
                                disp(conditionFields{condition})
                                disp(['Session: ', num2str(session)])
                                disp('   IN')
                                disp((squeeze(dataPoints_mean_in))')
                                disp('   NORM')
                                disp((squeeze(norm))')
                                disp('   OUT')
                                disp((squeeze(dataPoints_mean_out))')
                            end                           

                            
                            %normOut.(conditionFields{condition}).(dataTypes{field}).(auxParam{param})(session,:) = [];
                            %normOut.(conditionFields{condition}).(dataTypes{field}).(auxParam{param})(session,:) = []; % lower error bound
                            %normOut.(conditionFields{condition}).(dataTypes{field}).(auxParam{param})(session,:) = []; % upper error bound

                        end
                    end
                end                                    
            end
        end