function [dataOut, auxOut, auxOutPowers, subjects] = batch_pullOut_ERP(fileNameFields, outlierFilenameList, erpComponent, erpFilterType, handles)

    %% DEBUG
    debugMatFileName = 'tempPullOutERPs.mat';
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

    warning on
    
    numberOfFilesFound = length(fileNameFields);
    filesPerSubject = 12;
    numberOfSubjectsDone = numberOfFilesFound/filesPerSubject;
    disp([num2str(numberOfFilesFound), ' files found, thus most likely ', num2str(numberOfSubjectsDone), ' subjects done?'])
    
    if rem(numberOfSubjectsDone, 1) % if non-integer
       error('Lazy programmer have not managed yet the situation where you do not all the sessions done for each subject')       
    end
    
    % erpComponent
    outlierFilesFoundFromInput = false(length(outlierFilenameList), 1);
    
    for i = 1 : numberOfFilesFound                       
        
        % check intensity
        if fileNameFields{i}.intensity == 0
            intensity = 'dark';
        elseif fileNameFields{i}.intensity == 10
            intensity = 'dim';
        elseif fileNameFields{i}.intensity == 40
            intensity = 'bright';
        else
            intensity
            a = fileNameFields{i}.intensity
            error('Error in the intensity field!')
        end

        % check session
        if fileNameFields{i}.session == 1
            session = 'session1';
        elseif fileNameFields{i}.session == 2
            session = 'session2';
        elseif fileNameFields{i}.session == 3
            session = 'session3';
        elseif fileNameFields{i}.session == 4
            session = 'session4';
        else
            session
            error('Error in the session field!')
        end          
        
        % load the actual data
        disp(['    ... load the data from file: "', fileNameFields{i}.fileName, '"'])
        dataIn = load(fullfile(handles.path.matFilesOut, fileNameFields{i}.fileName));
        %dataIn.ERP_components
        %dataIn.ERP_components.(erpFilterType).target
                
        %% ERPs
        
            % check the contents of input                
                errorWithTheInput = 0;
                try
                    dataIn;
                catch err
                    % err 
                    disp(['      - ', fileNameFields{i}.fileName, ' contains no data, re-run the MAIN_ANALYSIS.m?'])
                    errorWithTheInput = 1;
                end
                
                if errorWithTheInput == 0
                    try
                        dataIn.ERP_components;
                    catch err
                        % err
                        disp(['      - ', fileNameFields{i}.fileName, ' contains no ERP component field, re-run the MAIN_ANALYSIS.m?'])
                        errorWithTheInput = 2;
                    end
                end

                if errorWithTheInput == 0
                    try
                        dataIn.ERP_components.(erpFilterType);
                    catch err
                        %err
                        disp(['      - ', fileNameFields{i}.fileName, ' has no data for the chosen filtering (', erpFilterType, '), re-run the MAIN_ANALYSIS.m?'])
                        errorWithTheInput = 3;
                    end
                end

                if errorWithTheInput == 0
                    try
                        dataIn.ERP_components.(erpFilterType).ERP;  
                    catch err
                        % err
                        disp(['      - ', fileNameFields{i}.fileName, ' has no ERP field for the chosen filtering (', erpFilterType, '), re-run the MAIN_ANALYSIS.m?'])
                        errorWithTheInput = 4;
                    end
                end
                
            % check if the input has been manually marked to be an outlier
            extraString = '_analyzed.mat';
            outlierListTrue = strcmp(strrep(fileNameFields{i}.fileName, extraString, '.bdf'), outlierFilenameList);
            thisFileMarkedAsOutlier = logical(sum(outlierListTrue));
            if thisFileMarkedAsOutlier == 1
                disp(['        - ', fileNameFields{i}.fileName, ' was marked to be an outlier'])
            end
            
                % check also that all the files marked as outliers are
                % found from the input files, i.e. if you have a typo in
                % your outlier files, it will be never found and corrected,
                % while you might think that it is been corrected
                outlierFilesFoundFromInput(outlierListTrue) = 1;

            % use more intuitive variable names
            if errorWithTheInput == 0 && thisFileMarkedAsOutlier == 0
                target = dataIn.ERP_components.(erpFilterType).ERP.target;
                distracter = dataIn.ERP_components.(erpFilterType).ERP.distr;
                standard = dataIn.ERP_components.(erpFilterType).ERP.std;                

                % Get structure of the input
                chNames = fieldnames(target);                
                componentNames = fieldnames(target.(chNames{1}){1});
                componentStatFields = fieldnames(target.(chNames{1}){1}.(componentNames{1}));                                
                
                % get number of epochs
                noOfEpochs.target = length(target.(chNames{1}));
                noOfEpochs.distr = length(distracter.(chNames{1}));
                noOfEpochs.std = length(standard.(chNames{1}));
                
            else
                try
                    target = batch_fillComponentFieldWithNaNs(chNames, componentNames, componentStatFields, noOfEpochs.target);
                    disp(['          - filling the target/distracter/standard ERPs with NaNs'])
                catch err
                    err
                    disp(['fileIndex = ', num2str(i)])
                    disp('does not work actually if your first file is faulty, fix later')
                end
                distracter = batch_fillComponentFieldWithNaNs(chNames, componentNames, componentStatFields, noOfEpochs.distr);
                standard = batch_fillComponentFieldWithNaNs(chNames, componentNames, componentStatFields, noOfEpochs.std);
            end
               
            % subject
            subject = fileNameFields{i}.subject;
            subjects{i} = subject;

            % assign to output, the desired component
            dataOut.(intensity).(session).target.component.(subject).(erpFilterType) = target;
            dataOut.(intensity).(session).distracter.component.(subject).(erpFilterType) = distracter;
            dataOut.(intensity).(session).standard.component.(subject).(erpFilterType) = standard;
        
        %% AUX variables
        analyzed_aux = dataIn.analyzed_aux;
        aux_fieldNames = fieldnames(analyzed_aux);
        
        indexFieldsFoundAmplit = 0;
        indexFieldsFoundPSD = 0;
        
        for ind = 1 : length(aux_fieldNames)
            
            subject = fileNameFields{i}.subject;
            
            if strcmp(aux_fieldNames{ind}, 'Amplit') || strcmp(aux_fieldNames{ind}, 'PSD')
                                
                powerFields = fieldnames(dataIn.analyzed_aux.(aux_fieldNames{ind}));
                
                if strcmp(aux_fieldNames{ind}, 'Amplit')
                    for j = 1 : length(powerFields)
                        indexFieldsFoundAmplit = indexFieldsFoundAmplit + 1;
                        auxOutPowers.(intensity).(session).(subject).Amplit.(powerFields{indexFieldsFoundAmplit}) = analyzed_aux.Amplit.(powerFields{indexFieldsFoundAmplit});
                        
                    end
                    
                elseif strcmp(aux_fieldNames{ind}, 'PSD')
                    
                    for j = 1 : length(powerFields)
                        indexFieldsFoundPSD = indexFieldsFoundPSD + 1;
                        auxOutPowers.(intensity).(session).(subject).PSD.(powerFields{indexFieldsFoundPSD}) = analyzed_aux.PSD.(powerFields{indexFieldsFoundPSD});
                    end
                    
                    % quick'n'dirty calculation of power ratios, ratio of
                    % alpha power compared to whole power as done by
                    % Cajochen et al. (2000) for example
                    % http://dx.doi.org/10.1016/S0166-4328(00)00236-9
                    auxOutPowers.(intensity).(session).(subject).ratio.alphaRatio = auxOutPowers.(intensity).(session).(subject).PSD.alpha / auxOutPowers.(intensity).(session).(subject).PSD.totalOzPz;
                    auxOutPowers.(intensity).(session).(subject).ratio.lowAlphaRatio = auxOutPowers.(intensity).(session).(subject).PSD.lowAlpha / auxOutPowers.(intensity).(session).(subject).PSD.totalOzPz;
                    auxOutPowers.(intensity).(session).(subject).ratio.highAlphaRatio = auxOutPowers.(intensity).(session).(subject).PSD.highAlpha / auxOutPowers.(intensity).(session).(subject).PSD.totalOzPz;
                    auxOutPowers.(intensity).(session).(subject).ratio.betaRatio = auxOutPowers.(intensity).(session).(subject).PSD.beta / auxOutPowers.(intensity).(session).(subject).PSD.totalFzCz;
                    auxOutPowers.(intensity).(session).(subject).ratio.deltaRatio = auxOutPowers.(intensity).(session).(subject).PSD.deltaThetaLockley / auxOutPowers.(intensity).(session).(subject).PSD.totalFzCz;
                    auxOutPowers.(intensity).(session).(subject).ratio.thetaRatio = auxOutPowers.(intensity).(session).(subject).PSD.theta / auxOutPowers.(intensity).(session).(subject).PSD.totalFzCz;

                    % According to Donskaya et al. (2012), ratio of theta /
                    % alpha shows a high correlation with subjective sleepiness
                    % than either alpha or theta power alone (end of 4st
                    % paragraph), http://dx.doi.org/10.1007/s11818-012-0561-1
                    auxOutPowers.(intensity).(session).(subject).ratio.alphaTheta = auxOutPowers.(intensity).(session).(subject).PSD.alpha / (auxOutPowers.(intensity).(session).(subject).PSD.alpha + auxOutPowers.(intensity).(session).(subject).PSD.theta);
                    auxOutPowers.(intensity).(session).(subject).ratio.lowAlphaTheta = auxOutPowers.(intensity).(session).(subject).PSD.lowAlpha / (auxOutPowers.(intensity).(session).(subject).PSD.lowAlpha + auxOutPowers.(intensity).(session).(subject).PSD.theta);
                    auxOutPowers.(intensity).(session).(subject).ratio.highAlphaTheta = auxOutPowers.(intensity).(session).(subject).PSD.highAlpha / (auxOutPowers.(intensity).(session).(subject).PSD.highAlpha + auxOutPowers.(intensity).(session).(subject).PSD.theta);
                    
                    % DEBUG
                    alpha = auxOutPowers.(intensity).(session).(subject).PSD.alpha;
                    beta = auxOutPowers.(intensity).(session).(subject).PSD.beta;
                    theta = auxOutPowers.(intensity).(session).(subject).PSD.theta;
                    delta = auxOutPowers.(intensity).(session).(subject).PSD.deltaThetaLockley;
                    FzCz_denom = beta + theta + delta;
                    total_FzCz = auxOutPowers.(intensity).(session).(subject).PSD.totalFzCz;
                    total_OzPz = auxOutPowers.(intensity).(session).(subject).PSD.totalOzPz;
                    
                    if (FzCz_denom / total_FzCz) > 1
                        warning('Ratios of beta+theta+delta should not exceed the total power in any conditions so you have a bug in the code probably?')
                    end
                    
                else
                    error('should not go here')
                end
            
            else                
                auxOut.(intensity).(session).(subject).scalars.(aux_fieldNames{ind}) = analyzed_aux.(aux_fieldNames{ind});
                
            end    
            
        end
     
    end % end of files
    
    %% check if all the outlier files have been found
    
        notFoundLinearIndices = find(outlierFilesFoundFromInput == 0);
        for i = 1 : length(notFoundLinearIndices)   
            warning('Not all the outlier files were found!')
            fprintf('     ')
            for j = 1 : i
               fprintf(' '); % add white space
            end        
            fprintf(['"', outlierFilenameList{notFoundLinearIndices(i)}, '" not found from the input files, is this correct?\n'])        
        end


