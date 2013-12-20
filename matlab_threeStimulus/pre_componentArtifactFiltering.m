function [matrixOut, matrixRegress, filtOnly, NaN_indices] = pre_componentArtifactFiltering(matrixIn, artifactIndices, rawMax, rowsIn, colsIn, loFreq, hiFreq, filterOrder, filterOrderSteep, dataType, parameters, handles)   

    
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempComponentArtifact.mat';
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
    

    %% Notch filter the mains (60 Hz) 
    % out using the modified filter from BioSig
    
        if strcmp(dataType, 'General') == 1
            matrixOut1a = zeros(rowsIn, colsIn);  % e.g. a lot of rows x 6 channels            
            HDR.SampleRate = parameters.EEG.srate;

                disp(['   remove 50/60 Hz (Notch, from BioSig package)'])
                for j = 1 : colsIn                        
                    matrixOut1a(:,j) = pre_remove5060hz_modified(matrixIn(:,j), HDR, 'PCA 60');                
                end
                dataIn = matrixOut1a; % for debug plotting
                filtOnly = matrixOut1a; 

        else
            matrixOut1a = matrixIn;        
            filtOnly = matrixOut1a; 
        end
        
        
    if parameters.artifacts.continuousFASTER == 1       
        
        % Bandpass filter
        [rowsIn, colsIn] = size(matrixOut1a);
        matrixOut1 = zeros(rowsIn, colsIn);  % e.g. a lot of rows x 6 channels            
        if strcmp(dataType, 'General')
            disp(['    - BANDPASS FILTERING'])
            disp(['        .. ', dataType,  ' (', num2str(loFreq), '-', num2str(hiFreq), ' Hz), order = ', num2str(filterOrder), ', ZeroPhase'])
            for j = 1 : colsIn                       
                matrixOut1(:,j) = pre_bandbassFilter(matrixOut1a(:,j), parameters.EEG.srate, [hiFreq, loFreq], filterOrder, filterOrderSteep, handles);                
            end
        end

        % Now split the continuous data into chunks and evaluate each chunk
        % as an epoch, and determine whether this chunk is artifacted or
        % not
        EEG = matrixOut1(:,1:handles.parameters.EEG.nrOfChannels);
        EOG = matrixOut1(:,handles.parameters.EEG.nrOfChannels+1);
        ECG = matrixOut1(:,handles.parameters.EEG.nrOfChannels+2);
        [matrixOut, matrixRegress, NaN_indices] = pre_FASTER_forContinuousData(EEG, EOG, ECG, parameters, handles);
        
        
        
    % "DUMB" OLD-SCHOOL FILTERING with EOG REGRESSION CORRECTION
    else
        
        
        % pre_continuousRecordingOldSchoolArtifactRemoval()
        
        %% Remove here fixed threshold artifacts (or get only the indices, 
        % and not pass the NaN values to the bandpass filter

            parameters.artifacts.continuousFASTER = 1;

            if strcmp(dataType, 'General') == 1
                [~, NaN_indices, numberOfNaNs] = pre_artifactFixedThreshold(matrixOut1a(:,1:handles.parameters.EEG.nrOfChannels), matrixOut1a(:,parameters.EEG.nrOfChannels + 1), parameters, handles);
                % noOfNans1 = sum(sum(isnan(matrixOut1a(:,1:handles.parameters.EEG.nrOfChannels)))) % there shouldn't be any at this point                
            else
                NaN_indices = artifactIndices;
            end

            %% Correct for EOG/ECG trends
            if strcmp(dataType, 'General') == 1

                EEGchans = parameters.EEG.nrOfChannels;            

                % You may want to consider later if you want to do the
                % regression for all types of filterings right here, and not
                % just for the GENERAL one that is going to be used for example
                % for TIME-FREQUENCY ANALYSIS, now the ERP method used for
                % assessing amplitude changes in ERP has not been pushed
                % through the regression correction

                % EOG
                if parameters.artifacts.applyRegressEOG == 1   
                    EOG = matrixOut1a(:,parameters.EEG.nrOfChannels + 1);
                    hiFreqEOG = parameters.artifacts.regressBirectContamCutoffs(2);
                    loFreqEOG = parameters.artifacts.regressBirectContamCutoffs(1);                
                    disp(['       EOG REGRESSION'])
                    disp(['        .. mitigate bidirectional contamination by filtering EOG (', num2str(loFreqEOG), '-', num2str(hiFreqEOG), ' Hz) before regression'])            
                    EOG = pre_bandbassFilter(EOG, parameters.EEG.srate, [hiFreqEOG, loFreqEOG], filterOrder, [], handles);   
                    EEG_withEOGcorr = pre_artifactRegressionWrapper(matrixOut1a, EOG, EEGchans, parameters);
                    disp(['          .. correct for EOG channel (regress_eog)'])
                end

                % ECG
                disp(['       ECG REGRESSION'])
                % not working very well atm at least (PT, 10 Dec 2013)
                % work on this later, see e.g. http://fsl.fmrib.ox.ac.uk/eeglab/fmribplugin/
                % and http://www.researchgate.net/publication/216221679_Matching_Pursuit_based_removal_of_cardiac_pulse-related_artifacts_in_EEGfMRI

                if parameters.artifacts.applyRegressECG == 1 && parameters.artifacts.applyRegressEOG == 1                
                    %EEG_withECGcorr = pre_artifactRegressionWrapper(EEG_withEOGcorr, matrixOut1a(:,parameters.EEG.nrOfChannels + 2), EEGchans, parameters);
                    disp(['        .. NOT CORRECTING ATM for ECG channel (pulse-related artifact)'])
                    EEG_withECGcorr = EEG_withEOGcorr;
                    matrixOut1a(:,1:parameters.EEG.nrOfChannels) = EEG_withECGcorr;


                elseif parameters.artifacts.applyRegressECG == 1                
                    %EEG_withECGcorr = pre_artifactRegressionWrapper(matrixOut1a, matrixOut1a(:,parameters.EEG.nrOfChannels + 2), EEGchans, parameters);
                    disp(['        .. NOT CORRECTING ATM for ECG channel (pulse-related artifact)'])
                    EEG_withECGcorr = EEG_withEOGcorr;
                    matrixOut1a(:,1:parameters.EEG.nrOfChannels) = EEG_withECGcorr;                

                else
                    disp(['        .. skip correction of ECG channel (pulse-related artifact)'])
                    matrixOut1a(:,1:parameters.EEG.nrOfChannels) = EEG_withEOGcorr;
                    EEG_withECGcorr = EEG_withEOGcorr;
                end

                % Output
                if parameters.artifacts.applyRegressECG == 1 || parameters.artifacts.applyRegressEOG == 1
                    matrixRegress = EEG_withECGcorr;
                else
                    matrixRegress = [];
                end                       

                %% PLOT the REGRESSION

                    t = linspace(1, length(matrixOut1a(:,1)), length(matrixOut1a(:,1))) / parameters.EEG.srate;
                    ch = 2; % 2 for Fz, 3 for Pz             
                    zoomRange = [160 165]; % in seconds
                    plot_artifactRegressionPlot(t, dataIn, EEG_withEOGcorr, EEG_withECGcorr, matrixRegress, ch, zoomRange, parameters, handles)           

            else
                filtOnly = matrixIn;
                matrixRegress = matrixIn;
            end

        %% Band-pass filter the data using a a subfunction            

            if strcmp(dataType, 'General')
                disp(['    - BANDPASS FILTERING'])
            end
            matrixOut1 = zeros(rowsIn, colsIn);  % e.g. a lot of rows x 6 channels            

            disp(['        .. ', dataType,  ' (', num2str(loFreq), '-', num2str(hiFreq), ' Hz), order = ', num2str(filterOrder), ', ZeroPhase'])
            for j = 1 : colsIn           
                if strcmp(dataType, 'General')
                    if parameters.artifacts.applyRegressECG == 1 || parameters.artifacts.applyRegressEOG == 1
                        matrixRegress = pre_bandbassFilter(matrixRegress, parameters.EEG.srate, [hiFreq, loFreq], filterOrder, filterOrderSteep, handles);                 
                        if j == 1
                            disp(['         .. filter the regression-corrected EEG as well'])
                        end
                    end
                end
                matrixOut1(:,j) = pre_bandbassFilter(matrixOut1a(:,j), parameters.EEG.srate, [hiFreq, loFreq], filterOrder, filterOrderSteep, handles);                
            end


        %% OUTPUT

            if strcmp(dataType, 'General') == 1
                % Don't convert the fixed threshold artifacts to NaNs
                matrixOut = matrixOut1;

            else
                % remove now the fixed threshold artifacts and do the other
                % operations
                % matrixOut(NaN_indices_corr) = NaN;
                % Don't convert the fixed threshold artifacts to NaNs
                matrixOut = matrixOut1;
            end
    end