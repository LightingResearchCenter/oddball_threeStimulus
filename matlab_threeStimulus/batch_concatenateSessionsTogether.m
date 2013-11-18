  function matricesConcatenated = batch_concatenateSessionsTogether(matricesIntensity, handles)

    ERPtypes = fieldnames(matricesIntensity);
    conditions = fieldnames(matricesIntensity.(ERPtypes{1}));
    noOfConditions = length(conditions);
    [noOfChannels, noOfSessions, noOfTrials, noOfSubjects] = size(matricesIntensity.(ERPtypes{1}).dark);

    matricesConcatenated = zeros(noOfChannels, noOfTrials*noOfSessions, noOfSubjects);

    for j = 1 : length(ERPtypes)
        for condition = 1 : noOfConditions
            for ch = 1 : noOfChannels
                for session = 1 : noOfSessions
                    for trials = 1 : noOfTrials                
                        sessionIndex = (session-1)*noOfTrials + trials;
                        for subjects = 1 : noOfSubjects
                            dataPoint = matricesIntensity.(ERPtypes{j}).(conditions{condition})(ch, session, trials, subjects);                  
                            matricesConcatenated.(ERPtypes{j}).(conditions{condition})(ch, sessionIndex, subjects) = dataPoint;
                            %{
                            plot(sessionIndex,dataPoint,'ko')
                            hold on
                            pause(0.2)
                            %}
                        end
                    end                
                end
            end
        end
    end

