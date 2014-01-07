function matricesIntensity = batch_getPerSubjectMatrix_Intensity(dataNorm, erpBandType, erpComponent, erpFilterType, fieldValue, noOfSubjects, noOfEpochs, handles)
       
    conditionFields = fieldnames(dataNorm);
    sessionFields = fieldnames(dataNorm.dark);
    noOfSessions = length(sessionFields);

    %{
    matricesIntensity.dark = zeros(handles.parameters.EEG.nrOfChannels, noOfSessions, noOfEpochs, noOfSubjects);
    matricesIntensity.dim = zeros(handles.parameters.EEG.nrOfChannels, noOfSessions, noOfEpochs, noOfSubjects);
    matricesIntensity.bright = zeros(handles.parameters.EEG.nrOfChannels, noOfSessions, noOfEpochs, noOfSubjects);
    %}

    dispNrOfCombinations = 0;
    for ii = 1 : length(conditionFields)                        

        for i = 1 : length(sessionFields)
            ERPtypes  = fieldnames(dataNorm.dark.(sessionFields{i}));

            for j = 1 : length(ERPtypes)
                subjectNames = fieldnames(dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType));

                for k = 1 : length(subjectNames)                    
                    try
                        chNames = fieldnames(dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType));
                    catch err
                        subjectNames{k}
                        fieldnames(dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}))
                        dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType)
                        fieldnames(dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType))
                        err
                        err.identifier                        
                    end

                    for l = 1 : length(chNames)
                        noOfEpochs = length(dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType).(chNames{l}));

                        for lj = 1 : noOfEpochs
                            try
                                compNames = fieldnames(dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType).(chNames{l}){lj});
                            catch err
                                [ii i j k l lj]
                                err
                                filterTypes = dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{1}).(erpBandType).(subjectNames{k})
                                chNames = dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{1}).(erpBandType).(subjectNames{k}).(erpFilterType)
                                epochs = dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{1}).(erpBandType).(subjectNames{k}).(erpFilterType).(chNames{l})
                            end

                            for m = 1 : 1 % length(compNames), we are only picking one

                                for n = 1 : 1 % length(statFields) we are only picking one                                        

                                    if dispNrOfCombinations == 0
                                        numberOfCombinations = length(conditionFields) * length(sessionFields) * length(ERPtypes) * length(subjectNames) * length(chNames) * noOfEpochs * 1 * 1;
                                        disp(['  .. .   Number of combinations: ', num2str(numberOfCombinations)])
                                        dispNrOfCombinations = 1;
                                    end                                        

                                    % easier variable names
                                    if ~strcmp(erpComponent, 'RT')
                                        % disp([i lj k])
                                        try
                                            dataPoint = dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType).(chNames{l}){lj}.(erpComponent).(fieldValue);
                                        catch err
                                            err
                                            erpComponent
                                            subjectNames{k}
                                            dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType)
                                            dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType).(chNames{l}){lj}
                                            ataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType).(chNames{l}){lj}.(erpComponent)
                                            error('Typo with what you wanted to get?')
                                        end
                                        if isstruct(dataPoint)
                                            error('Data point at this part should not be a structure, error in the code! Check how structures are nested and maybe you forgot one field?')
                                        end
                                        matricesIntensity.(ERPtypes{j}).(conditionFields{ii})(l,i,lj,k) = dataPoint;
                                    else      
                                        try 
                                            dataPoint = dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType).(chNames{l}){lj}.(erpComponent);
                                        catch err
                                            err
                                            a1 = dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType)
                                            a2 = dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType).(chNames{l})
                                            whos
                                            a3 = dataNorm.(conditionFields{ii}).(sessionFields{i}).(ERPtypes{j}).(erpBandType).(subjectNames{k}).(erpFilterType).(chNames{l}){lj}
                                            whos
                                        end
                                        if isstruct(dataPoint)
                                            error('Data point at this part should not be a structure, error in the code! Check how structures are nested and maybe you forgot one field? Maybe an extra fieldValue for RT from previous steps?')
                                        end
                                        matricesIntensity.(ERPtypes{j}).(conditionFields{ii})(l,i,lj,k) = dataPoint;
                                            

                                    end
                                end
                            end
                        end                                
                    end                        
                end
            end
        end
    end
