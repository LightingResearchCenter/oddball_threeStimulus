 function matricesStat = batch_averageTheSubjectsMatrix(matricesConcatenated, statsPer, handles)
        
    ERPtypes = fieldnames(matricesConcatenated);
    conditions = fieldnames(matricesConcatenated.(ERPtypes{1}));
    noOfConditions = length(conditions);
    [noOfChannels, noOfTrials, noOfSubjects] = size(matricesConcatenated.(ERPtypes{1}).dark);

    dim = 2;
    flag = 0;

    for j = 1 : length(ERPtypes)
        for condition = 1 : noOfConditions
            for ch = 1 : noOfChannels                
                for subjects = 1 : noOfSubjects
                    if strcmp(statsPer, 'session') 
                        matricesStat.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}) = batch_calculateStatsForMatrix(matricesConcatenated.(conditions{condition})(ch, :, :), dim, flag, handles);
                    elseif strcmp(statsPer, 'trials') 
                        matricesStat.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}) = batch_calculateStatsForMatrix(matricesConcatenated.(conditions{condition})(ch, :, :), dim, flag, handles);
                    end
                end            
            end
        end
    end