function matrixNorm = batch_normalizeComponents(matricesSessionAver, typeOfNorm, subjects, handles)
    
    %% DEBUG
    debugMatFileName = 'tempAverFirstSession.mat';
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
        normValues = batch_normalizeSubLoop(matricesSessionAver, [], typeOfNorm, [], 'dark', 'getNorm', subjects, handles);
        matrixNorm = batch_normalizeSubLoop(matricesSessionAver, normValues, typeOfNorm, [], 'dark', 'normalize', subjects, handles);
    elseif strcmp(typeOfNorm, 'firstSession')
        % Go through the data
        normValues = batch_normalizeSubLoop(matricesSessionAver, [], typeOfNorm, 1, [], 'getNorm', subjects, handles);
        matrixNorm = batch_normalizeSubLoop(matricesSessionAver, normValues, typeOfNorm, 1, [], 'normalize', subjects, handles);
    end

    
    
    function normOut = batch_normalizeSubLoop(matricesSessionAver, normValues, typeOfNorm, sessionToFix, conditionToFix, normPassType, subjects, handles)
        
        normOut = matricesSessionAver;
        
        % Get field parameters
        ERPtypes = fieldnames(matricesSessionAver);
        conditions = fieldnames(matricesSessionAver.(ERPtypes{1}));
        noOfConditions = length(conditions);
        noOfChannels = length(fieldnames(matricesSessionAver.(ERPtypes{1}).(conditions{1})));        
        noOfTrialsOrSessions = length(matricesSessionAver.(ERPtypes{1}).(conditions{1}).Fz);
        statFields = fieldnames(matricesSessionAver.(ERPtypes{1}).(conditions{1}).Fz{1});
        noOfSubjects = length(matricesSessionAver.(ERPtypes{1}).(conditions{1}).Fz{1}.mean);   
        
        if strcmp(typeOfNorm, 'darkCondition')
            conditionFixedIndex = find(strcmp(conditions, conditionToFix) == 1);
        elseif strcmp(typeOfNorm, 'firstSession')
            sessionFixedIndex = sessionToFix;
        else
            error('typo?')
        end
        
        for condition = 1 : noOfConditions        
            for j = 1 : length(ERPtypes)
                for ch = 1 : noOfChannels
                    for sessionTrial = 1 : noOfTrialsOrSessions
                        % for subjects = 1 : noOfSubjects
                            % for stat = 1 : length(statFields)
                                
                                if strcmp(normPassType, 'getNorm')                                    
                                    
                                    if strcmp(typeOfNorm, 'darkCondition') 
                                        % keep the index fixed for the
                                        % condition that you save to the output
                                        % cell structure
                                        dataPoints = matricesSessionAver.(ERPtypes{j}).(conditions{conditionFixedIndex}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){sessionTrial}.mean;                                        
                                        normOut.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){sessionTrial}.mean = dataPoints;
                                        
                                        if ch == 1 && strcmp(ERPtypes{j}, 'target') && 1 == 2
                                            disp(conditions{condition})
                                            disp([j ch sessionTrial dataPoints])
                                        end
                                        
                                    elseif strcmp(typeOfNorm, 'firstSession')                                        
                                        dataPoints = matricesSessionAver.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){sessionFixedIndex}.mean;
                                        normOut.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){sessionTrial}.mean = dataPoints;
                                    end   

                                elseif strcmp(normPassType, 'normalize')
                                    
                                    % does not matter know on what you have
                                    % decided to normalize the data to,
                                    % just matters that the first pass was
                                    % correct
                                                                
                                    dataPoints_mean_in = matricesSessionAver.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){sessionTrial}.mean;
                                    dataPoints_SD_in = matricesSessionAver.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){sessionTrial}.SD;
                                    norm = normValues.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){sessionTrial}.mean;
                                    
                                    [dataPoints_mean_out, ~, ~] = batch_normalizeLowLevel(dataPoints_mean_in, 'mean', [], [], norm, handles);                                         
                                    [dataPoints_SD_out, LE, UE] = batch_normalizeLowLevel(dataPoints_SD_in, 'SD', dataPoints_mean_in, dataPoints_mean_out, norm, handles);                                    
                                                                        
                                    % exclude outliers
                                    [dataPoints_mean_out, outlierIndices] = batch_excludeOutliersDuringBatch(dataPoints_mean_out, subjects, 'normalize', handles);
                                    
                                    % assign out
                                    normOut.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){sessionTrial}.mean = dataPoints_mean_out;
                                    
                                    % use the outlier indices for other
                                    % variables
                                    dataPoints_SD_out(outlierIndices) = NaN;
                                    LE(outlierIndices) = NaN;
                                    UE(outlierIndices) = NaN;
                                                                        
                                    % assign out
                                    normOut.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){sessionTrial}.SD = dataPoints_SD_out;
                                    normOut.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){sessionTrial}.LE = LE; % lower error bound
                                    normOut.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}){sessionTrial}.UE = UE; % upper error bound
                                    
                                    if sum(~isnan(dataPoints_mean_out)) > 0 && 1 == 2
                                        disp('             normalize')    
                                        disp(dataPoints_mean_in)
                                        disp(dataPoints_SD_in)
                                        disp(norm)
                                        disp(dataPoints_mean_out)
                                        disp(dataPoints_SD_out)
                                        
                                    end
                                    
                                    if ch == 1 && strcmp(ERPtypes{j}, 'target') && 1 == 2                                        
                                        disp('             normalize') 
                                        disp(conditions{condition})
                                        disp([j ch sessionTrial dataPoints_mean_in])
                                    end
                                    
                                end 
                            % end                            
                        % end
                    end
                end
            end
        end

        
    %% AVERAGE OVER SUBJECTS
    
    % different ways to normalize the data

        % note that the dimension of the data should not be changed
        % during the normalization (only some normalization value are
        % obtained and then one divides the input with that) so the
        % different "primitive normalizations" can be combined if
        % wanted

        % 1) All the dark values are zero

            % Average over subjects so that each individual is first normalized
            % based on the dark session
            %{
            disp('  .. normalize for dark condition')
            normData = batch_normalizeForDarkCondition(dataOut, 'getNorm', [], erpComponent, erpDataType, fieldValue, handles);
            dataNorm.darkCondition = batch_normalizeForDarkCondition(dataOut, 'normalize', normData, erpComponent, erpDataType, fieldValue, handles);

        % 2) Normalize to first trial

            % The first session (pre-light baseline condition) is
            % normalized to unity (just change the index, and you can
            % normalize to last session, or 2nd, or 3rd, etc.)
            disp('   .. normalize for first session')
            sessionIndexForNorm = 1;
            normData = batch_normalizeForGivenSession(dataOut, 'getNorm', sessionIndexForNorm, [], erpComponent, erpDataType, fieldValue, handles);
            dataNorm.firstSession = batch_normalizeForGivenSession(dataOut, 'normalize', [], normData, erpComponent, erpDataType, fieldValue, handles);
            %}

        % 3) From the first trial, average for the dark session

            %{
            % see what is the bug
            disp('    .. normalize for dark condition from "first session" normalization')
            normData = batch_normalizeForDarkCondition(dataNorm.firstSession, 'getNorm', [], erpComponent, erpDataType, fieldValue, handles);
            dataNorm.darkAfterFirst = batch_normalizeForDarkCondition(dataOut, 'normalize', normData, erpComponent, erpDataType, fieldValue, handles);   
            %}

        % 4) RAW Without any normalization
            % disp('     .. raw without any normalization')
            % dataNorm.rawWithNoNormalization = dataOut;
                