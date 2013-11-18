function dataOut = batch_extraCreateMatrix(dataIn, sensorType, scalarOrVector, handles)

    conditionFields = fieldnames(dataIn);
    sessionFields = fieldnames(dataIn.(conditionFields{1}));
    subjectFields = fieldnames(dataIn.(conditionFields{1}).(sessionFields{1}));
    scalarFields = fieldnames(dataIn.(conditionFields{1}).(sessionFields{1}).(subjectFields{1}).scalar);
    vectorFields = fieldnames(dataIn.(conditionFields{1}).(sessionFields{1}).(subjectFields{1}).vector);

    noOfSessions = length(sessionFields);
    noOfSubjects = length(subjectFields);

    if strcmp(scalarOrVector, 'vector')
        numberOfFields = length(vectorFields);
        fields = vectorFields;        
    elseif strcmp(scalarOrVector, 'scalar')
        numberOfFields = length(scalarFields);
        fields = scalarFields;
    else
        error('typo here probably, variable scalarOrVector should be either ''scalar'' or ''vector''')
    end
      

    for intensity = 1 : length(conditionFields)
        for session = 1 : noOfSessions
            for subject = 1 : noOfSubjects
                for field = 1 : numberOfFields    
                    
                    % for EEG data, this is cell, each channel as cell
                    % element
                    dataIn.(conditionFields{intensity}).(sessionFields{session}).(subjectFields{subject}).(scalarOrVector).(fields{field});
                    isCellOrNot = iscell(dataIn.(conditionFields{intensity}).(sessionFields{session}).(subjectFields{subject}).(scalarOrVector).(fields{field}));
                    
                    if isCellOrNot == 1
                        
                        for ch = 1 : length(dataIn.(conditionFields{intensity}).(sessionFields{session}).(subjectFields{subject}).(scalarOrVector).(fields{field}))
                           
                            dataOut.(conditionFields{intensity}).(scalarOrVector).(fields{field}){ch}(:, session, subject) = ...
                                dataIn.(conditionFields{intensity}).(sessionFields{session}).(subjectFields{subject}).(scalarOrVector).(fields{field}){ch};
                            
                        end

                    else
                        %{
                        dataIn
                        dataIn.(conditionFields{intensity})
                        dataIn.(conditionFields{intensity}).(sessionFields{session})
                        dataIn.(conditionFields{intensity}).(sessionFields{session}).(subjectFields{subject})
                        dataIn.(conditionFields{intensity}).(sessionFields{session}).(subjectFields{subject}).(scalarOrVector)

                        a = fields{field}
                        %}

                        dataOut.(conditionFields{intensity}).(scalarOrVector).(fields{field})(:, session, subject) = ...
                            dataIn.(conditionFields{intensity}).(sessionFields{session}).(subjectFields{subject}).(scalarOrVector).(fields{field});
                    end                      

                    
                end
            end
        end
    end
    