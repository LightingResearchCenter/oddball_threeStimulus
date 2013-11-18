function epochsOut = pre_removeBaseline(epochs, parameters, handles)

    baselineCorr = (parameters.oddballTask.ERP_baseline * parameters.EEG.srate);

    for i = 1 : length(epochs.oddball_regular)       
        epochsOut.oddball_regular{i} = epochs.oddball_regular{i}(baselineCorr+1:end,:);
    end        
        
    for i = 1 : length(epochs.oddball_irregular)
        epochsOut.oddball_irregular{i} = epochs.oddball_irregular{i}(baselineCorr+1:end,:);
    end

    