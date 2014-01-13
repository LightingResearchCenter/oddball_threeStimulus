function statsTests = batch_calculateStatSignificances(statsOut, normFieldName, erpBandType, erpComponent, erpFilterType, matrixIn, subjects, dataType, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempStatSignif.mat';
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
    end

    % example input data
    if strcmp(dataType, 'component')
        noOfSessions = size(statsOut.target.dark.Fz.mean);    
        [noOfSessions, noOfSubjects] = size(statsOut.target.dark.Fz.meanIn);
        statFields = fieldnames(statsOut.target.dark.Fz)
    elseif strcmp(dataType, 'AUX')
        noOfSessions = length(statsOut.dark.PSD.alpha.mean);
        [noOfSessions, noOfSubjects] = size(matrixIn.dark.PSD.alpha);
        %.PSD.alpha.mean
    elseif strcmp(dataType, 'extraHeart')
        matrixIn = squeeze(matrixIn); % why is there the extra singleton dimension actually?
        [noOfSessions, noOfSubjects] = size(matrixIn.dark.scalar.HR_Mean);
    else
        error(['Your datatype is: ', dataType, '. There is nothing coded for this dataType, or is this a typo?'])
        
        
    end


    % construct the comparisons
    
    if strcmp(dataType, 'component')
        stimTypes = fieldnames(statsOut);
        conditions = fieldnames(statsOut.(stimTypes{1}));
        chNames = fieldnames(statsOut.(stimTypes{1}).(conditions{1}));

        statField = 'meanIn';

        for stim = 1 : length(stimTypes)

            for ch = 1 : length(chNames)

                for cond = 1 : length(conditions)
                    comparisonMatrix(cond,:,:) = statsOut.(stimTypes{stim}).(conditions{cond}).(chNames{ch}).(statField)                
                end

                [statsTests.assumptions.(stimTypes{stim}).(chNames{ch}), statsTests.testResults.(stimTypes{stim}).(chNames{ch})] = ...
                    stat_signifTestWrapper(comparisonMatrix, statsOut.(stimTypes{stim}), subjects, normFieldName, erpBandType, erpComponent, erpFilterType, statField, stimTypes{stim}, chNames{ch}, dataType, handles.parameters, handles)

            end
        end
        
    elseif strcmp(dataType, 'AUX')
        
        conditions = fieldnames(statsOut);
        statField = 'mean';
        
        powerTypes = fieldnames(statsOut.dark); % PSD / amplit / ratio
         
        for powerType = 1 : length(powerTypes)
            
            bandRatioTypes = fieldnames(statsOut.dark.(powerTypes{powerType})); % alpha, beta, etc.
            for bandType = 1 : length(powerTypes)
               
                for cond = 1 : length(conditions)            
                    matrixInputPerCondition = matrixIn.(conditions{cond}).(powerTypes{powerType}).(bandRatioTypes{bandType});
                    comparisonMatrix(cond,:,:) = matrixInputPerCondition;
                end
                
                [statsTests.assumptions.(powerTypes{powerType}).(bandRatioTypes{bandType}), statsTests.testResults.(powerTypes{powerType}).(bandRatioTypes{bandType})] = ...
                    stat_signifTestWrapper(comparisonMatrix, statsOut, subjects, normFieldName, erpBandType, erpComponent, erpFilterType, statField, powerTypes{powerType}, bandRatioTypes{bandType}, dataType, handles.parameters, handles);
                
            end
            
        end
        
    elseif strcmp(dataType, 'extraHeart')
        
        conditions = fieldnames(statsOut);
        statField = 'mean';
        
        dimType = 'scalar';
        measureFields = fieldnames(matrixIn.(conditions{1}).(dimType));
        
        for meas = 1 : length(measureFields)
            
            for cond = 1 : length(conditions)            
                matrixInputPerCond = matrixIn.(conditions{1}).(dimType).(measureFields{meas});
                comparisonMatrix(cond,:,:) = matrixInputPerCond;
            end            
            
             [statsTests.assumptions.(measureFields{meas}), statsTests.testResults.(measureFields{meas})] = ...
                    stat_signifTestWrapper(comparisonMatrix, statsOut, subjects, normFieldName, [], [], [], statField, measureFields{meas}, [], dataType, handles.parameters, handles);
        
            
        end
        
       
        
    else
        
    end


