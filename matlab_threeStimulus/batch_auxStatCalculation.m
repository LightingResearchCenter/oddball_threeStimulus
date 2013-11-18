function auxStat = batch_auxStatCalculation(auxMatrix, handles)
        
    normFields = fieldnames(auxMatrix);
    dim = 2;
    flag = 0;
    
    for i = 1 : length(normFields)
        conditions = fieldnames(auxMatrix.(normFields{i}));
    
        for cond = 1 : length(conditions)            
            auxFields = fieldnames(auxMatrix.(normFields{i}).(conditions{cond}));
        
            for field = 1 : length(auxFields)
                auxParam = fieldnames(auxMatrix.(normFields{i}).(conditions{cond}).(auxFields{field}));
            
                for param = 1 : length(auxParam)
                    
                    [noOfSessions, noOfSubjects] = size(auxMatrix.(normFields{i}).(conditions{cond}).(auxFields{field}).(auxParam{param}));
                  
                    matrixIn = auxMatrix.(normFields{i}).(conditions{cond}).(auxFields{field}).(auxParam{param})(:,:);
                    structOut = batch_calculateStatsForMatrix(matrixIn, dim, flag, handles);
                    auxStat.(normFields{i}).(conditions{cond}).(auxFields{field}).(auxParam{param}) = structOut;

                    if 1 == 1 % debug

                        matrixIn
                        meanDebug = structOut.mean
                        sessionsMean = auxStat.(normFields{i}).(conditions{cond}).(auxFields{field}).(auxParam{param}).mean
                        sessionsSD = auxStat.(normFields{i}).(conditions{cond}).(auxFields{field}).(auxParam{param}).SD
                        pause
                    
                    end
                        
                    
                    
                end
            end
        end
    end

