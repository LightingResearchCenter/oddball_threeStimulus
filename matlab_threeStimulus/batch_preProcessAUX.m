function [auxStat, auxStatPowers] = batch_preProcessAUX(auxOut, auxOutPowers, handles)

    %% DEBUG
    debugMatFileName = 'tempPreProcessAuxScalars.mat';
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
    
    %% AUX OUT : Scalars like IAF, or whatever you think of putting inside
    
        matOut_aux = batch_auxCellToMatrix(auxOut, handles);
        
        % normalize 
        auxMatrix.absolute = matOut_aux;
        auxMatrix.darkCondition = batch_normalizeAux(matOut_aux, 'darkCondition', handles);
        auxMatrix.firstSession = batch_normalizeAux(matOut_aux, 'firstSession', handles);
        
        % calculate the stats subjects
        auxStat = batch_auxStatCalculation(auxMatrix, handles);        
    

    %% AUX OUT / POWERS : EEG Band Powers / Ratios
    
        matOut_auxPowers = batch_auxCellToMatrix(auxOutPowers, handles);        
        
        % normalize 
        auxMatPowers.absolute = matOut_auxPowers;
        auxMatPowers.darkCondition = batch_normalizeAux(matOut_auxPowers, 'darkCondition', handles);
        auxMatPowers.firstSession = batch_normalizeAux(matOut_auxPowers, 'firstSession', handles);
        
        % calculate the stats subjects
        auxStatPowers = batch_auxStatCalculation(auxMatPowers, handles);
    
    
    function matOut_aux = batch_auxCellToMatrix(auxOut, handles)

        conditions = fieldnames(auxOut);
        sessions = fieldnames(auxOut.(conditions{1}));
        subject = fieldnames(auxOut.(conditions{1}).(sessions{1}));
        auxFields = fieldnames(auxOut.(conditions{1}).(sessions{1}).(subject{1}));

        for i = 1 : length(conditions)
            for f = 1 : length(auxFields)
                auxParam = fieldnames(auxOut.(conditions{i}).(sessions{1}).(subject{1}).(auxFields{f}));
                for p = 1 : length(auxParam)
                    
                    for session = 1 : length(sessions)
                        for sub = 1 : length(subject)
                            
                            % disp([auxFields{f} ' ' auxParam{p}])
                            subjectDataPoint = auxOut.(conditions{i}).(sessions{session}).(subject{sub}).(auxFields{f}).(auxParam{p});
                            matOut_aux.(conditions{i}).(auxFields{f}).(auxParam{p})(session, sub) = subjectDataPoint;

                        end
                    end
                    
                end
            end
        end        

        % debug = matOut_aux.(conditions{i}).(auxFields{f}).(auxParam{p})
