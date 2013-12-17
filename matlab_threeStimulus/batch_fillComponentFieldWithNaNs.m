function erpType = batch_fillComponentFieldWithNaNs(chNames, componentNames, componentStatFields, noOfEpochs)

    for ch = 1 : length(chNames)
        for ep = 1 : noOfEpochs
            for comp = 1 : length(componentNames)
                for statField = 1 : length(componentStatFields)
                    erpType.(chNames{ch}){ep}.(componentNames{comp}).(componentStatFields{statField}) = NaN;                
                end            
            end        
        end
    end