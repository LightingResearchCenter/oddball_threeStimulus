function [dataOut, auxOut, auxOutPowers] = batch_pullOut_ERP(fileNameFields, erpComponent, erpDataType, handles)

    %% DEBUG
    debugMatFileName = 'tempDataOut.mat';
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

    numberOfFilesFound = length(fileNameFields);
    filesPerSubject = 12;
    numberOfSubjectsDone = numberOfFilesFound/filesPerSubject;
    disp([num2str(numberOfFilesFound), ' files found, thus most likely ', num2str(numberOfSubjectsDone), ' subjects done?'])
    
    if rem(numberOfSubjectsDone, 1) % if non-integer
       error('Lazy programmer have not managed yet the situation where you do not all the sessions done for each subject')       
    end
    
    % erpComponent
    
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
        %dataIn.ERPsOut
        %dataIn.ERPsOut.target
                
        %% ERPs
        
            % use more intuitive varriable name
            ERP_perChannels.target = dataIn.ERPsOut.target;
            ERP_perChannels.distracter = dataIn.ERPsOut.distracter;
            ERP_perChannels.standard = dataIn.ERPsOut.standard;

                % correct later the switch        

            % determine the desired data
            if strcmp(erpDataType, 'EP')            
                target = ERP_perChannels.target;
                distracter = ERP_perChannels.distracter;
                standard = ERP_perChannels.standard;

            elseif strcmp(erpDataType, 'filt')
                target = ERP_perChannels.target;
                distracter = ERP_perChannels.distracter;
                standard = ERP_perChannels.standard;

            elseif strcmp(erpDataType, 'CNV')
                warning('NOT IMPLEMENTED YET, and NOT SAVED from PROCESS_singleFile')

            elseif strcmp(erpDataType, 'Alpha')
                warning('NOT IMPLEMENTED YET, and NOT SAVED from PROCESS_singleFile')

            elseif strcmp(erpDataType, 'General')
                warning('NOT IMPLEMENTED YET, and NOT SAVED from PROCESS_singleFile')

            else
                erpDataType
                error(['erpDataType = ', erpDataType, 'What kind of ERP data type you actually wanted to process, maybe a typo?'])
            end

            subject = fileNameFields{i}.subject;

            % assign to output, the desired component
            dataOut.(intensity).(session).target.component.(subject) = target;
            dataOut.(intensity).(session).distracter.component.(subject) = distracter;
            dataOut.(intensity).(session).standard.component.(subject) = standard;
        
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
     
    end
    
    