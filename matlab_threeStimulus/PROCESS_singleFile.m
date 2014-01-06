function [epochs, analyzed, TF, dataMatrix_filtGeneral, alpha, powers, handles] = ...
                PROCESS_singleFile(inputFiles, fileNameIn, dataMatrixIn, triggers, handles)

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempProcess.mat';
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
    
    disp(' ')                
    analyzed = [];
    TF = [];

    %% PRE-PROCESS THE DATA            
    [dataMatrix_filtGeneral, dataMatrix_filtGeneralContinuous, firstBandpassFilteredMatrix, artifactNaN_indices, ...
        dataMatrix_filtAlpha, dataMatrix_filt, dataMatrix_filtSmooth1, dataMatrix_filtSmooth2, dataMatrix_filt_CNV, dataMatrix_filtP300] = ...
        process_preProcessFiltering(dataMatrixIn, handles);
    
            
    %% GENERAL Time-Series Analysis
    disp(' ')
    disp('    Processing Time-series EEG')
        
        EEGind = [1 handles.parameters.EEG.nrOfChannels];
        
        % raw still contain the reference channels
        EOGind = [handles.parameters.EEG.nrOfChannels+1 handles.parameters.EEG.nrOfChannels+1];
        EOGind_raw = [handles.parameters.EEG.nrOfChannels+1+handles.parameters.BioSemi.chOffset handles.parameters.EEG.nrOfChannels+1+handles.parameters.BioSemi.chOffset];
        
        % raw still contain the reference channels
        ECGind = [handles.parameters.EEG.nrOfChannels+2 handles.parameters.EEG.nrOfChannels+2];
        ECGind_raw = [handles.parameters.EEG.nrOfChannels+2+handles.parameters.BioSemi.chOffset handles.parameters.EEG.nrOfChannels+2+handles.parameters.BioSemi.chOffset];

        % PROCESS the "Time-Series EEG", i.e. al the channels without
        % ERP oddball epoching, artifacts removed and bandpass-filtered     
        handles.parameters.compute_MFDFA = 0;
        [alpha, powers, amplitSpectrum, PSD, SEM, heart, fractalAnalysis] =  process_timeSeriesEEG(dataMatrix_filtGeneralContinuous(:,EEGind(1):EEGind(2)), ... % EEG
                                                                            dataMatrix_filtGeneralContinuous(:,EOGind(1):EOGind(2)), ... % EOG 
                                                                            dataMatrixIn(:,EOGind_raw(1):EOGind_raw(2)), ... % EOG RAW
                                                                            dataMatrix_filtGeneralContinuous(:,ECGind(1):ECGind(2)), ... % ECG
                                                                            dataMatrixIn(:,ECGind_raw(1):ECGind_raw(2)), ... % ECG RAW
                                                                            triggers, handles.style, handles.parameters, handles);                                                                           
        close all
                       
    %% Epoch the EEG
    
        % i.e. Split into Oddball, Distracter and Standard Event-Related Potentials (ERPs)            
        disp('    Split recording to ERP epochs')              
        
        filteringType = 'bandpass';        

            ERP_baseline = handles.parameters.oddballTask.ERP_baseline;
            ERP_duration = handles.parameters.oddballTask.ERP_duration;
            
            % Now all the different possible filtering/data types are
            % stored to fields of a structure, and you can rather freely
            % add new types and then just use that added type easily from
            % the BATCH portion of the code without rewriting any code
        
            erpType = 'ERP'; disp(['     ', erpType]) % use the regression corrected
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices] = ...
                pre_epochToERPs(dataMatrix_filt(:,:), triggers, [], ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles);    

            erpType = 'RAW'; disp(['     ', erpType])
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices] = ...
                pre_epochToERPs(dataMatrixIn(:,:), triggers, epoch_indices, ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles); 
            
            erpType = 'ERPsmooth1'; disp(['     ', erpType])
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices] = ...
                pre_epochToERPs(dataMatrix_filtSmooth1(:,:), triggers, epoch_indices, ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles);                              
            
            erpType = 'ERPsmooth2'; disp(['     ', erpType])
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices] = ...
                pre_epochToERPs(dataMatrix_filtSmooth2(:,:), triggers, epoch_indices, ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles);                              
            
            erpType = 'GENERAL'; disp(['     ', erpType])
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices] = ...
                pre_epochToERPs(dataMatrix_filtGeneral(:,:), triggers, epoch_indices, ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles);
            
            erpType = 'GENERAL_regress'; disp(['     ', erpType])
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices] = ...
                pre_epochToERPs(dataMatrix_filtGeneralRegress(:,:), triggers, epoch_indices, ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles);
                        
            erpType = 'CNV'; disp(['     ', erpType])
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices] = ...                
                pre_epochToERPs(dataMatrix_filt_CNV(:,:), triggers, epoch_indices, ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles);
                        
            erpType = 'ALPHA'; disp(['     ', erpType])
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices] = ...
                pre_epochToERPs(dataMatrix_filtAlpha(:,:), triggers, epoch_indices, ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles);
                        
            erpType = 'P300'; disp(['     ', erpType])
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices] = ...
                pre_epochToERPs(dataMatrix_filtP300(:,:), triggers, epoch_indices, ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles);
            
            erpType = 'ECG'; disp(['     ', erpType])
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices] = ...
                pre_epochToERPs(dataMatrixIn(:,handles.parameters.EEG.nrOfChannels+2:handles.parameters.EEG.nrOfChannels+2), triggers, epoch_indices, ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles);
                        
            erpType = 'EOG'; disp(['     ', erpType])
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices] = ...
                pre_epochToERPs(dataMatrix_filtGeneral(:,handles.parameters.EEG.nrOfChannels+1:handles.parameters.EEG.nrOfChannels+1), triggers, epoch_indices, ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles);
            
            % Note! Now the ERP for EP Denoising need to have equal lengths
            % for pre-stimulus baseline and post-stimulus recording which
            % is not going to be case (necessarily) for the other epochs
            ERP_baseline = handles.parameters.oddballTask.ERP_denoisingEP_epochLimits(1) * -1;
            ERP_duration = handles.parameters.oddballTask.ERP_denoisingEP_epochLimits(2);
            
            % Note that now, the epochs only match the EP_DEN requirements,
            % and data is still the same as with ERP bandpass-filtered, and
            % we are going to do the actual denoising below
            erpType = 'EP_DEN'; disp(['     ', erpType])
            [epochs.(filteringType).(erpType).target, epochs.(filteringType).(erpType).distr, epochs.(filteringType).(erpType).std, epoch_indices_EP] = ...
                pre_epochToERPs(dataMatrix_filt(:,:), triggers, [], ERP_baseline, ERP_duration, alpha, handles.parameters, erpType, handles);
            
            % release some memory
            clear dataMatrix    
            clear dataMatrix_filt
            clear ref
            clear triggersRaw
            clear triggers                 

            % Save into GDF?
            % http://en.wikipedia.org/wiki/General_Data_Format_for_Biomedical_Signals                       
    
    %% FASTER: Artifact removal for the concatenad RAW / GENERALLY filtered

        if handles.parameters.artifacts.useFASTER == 1
            
            % use the ERP-filtered epochs as the reference as it should be
            % smoother      
            disp(' ')
            disp('    FASTER Artifact Rejection')      
            
            filteringType = 'bandpass';
            rejectBasedOn = 'GENERAL'; % FASTER paper recommends some filtering before the algorithm
            referenceType = 'RAW'; % for plotting etc.
            debugOn = 1; % show debug plots
            
            % Use a subfunction wrapper (with a lot of inputs and outputs)
            [epochs.(filteringType).(rejectBasedOn).target, epochs.(filteringType).(rejectBasedOn).distr, epochs.(filteringType).(rejectBasedOn).std,...
                artifactIndices_FASTER.target, artifactIndices_FASTER.distr, artifactIndices_FASTER.std]...
                = pre_FASTER_forAllStimuli(epochs.(filteringType).(rejectBasedOn).target, epochs.(filteringType).(rejectBasedOn).distr, epochs.(filteringType).(rejectBasedOn).std, ... % reject
                                           epochs.(filteringType).(referenceType).target, epochs.(filteringType).(referenceType).distr, epochs.(filteringType).(referenceType).std, ... % reference
                                           epochs.(filteringType).EOG.target, epochs.(filteringType).EOG.distr, epochs.(filteringType).EOG.std, ... % EOG
                                           epochs.(filteringType).ECG.target, epochs.(filteringType).ECG.distr, epochs.(filteringType).ECG.std, ... % ECG
                                           filteringType, rejectBasedOn, referenceType, ... % strings 
                                           debugOn, ... % boolean flags
                                           handles.parameters, handles); % parameters                                       
                                       
            % Now we have found the artifact indices based on the RAW
            % (unfiltered epochs) and we need to use the found indices to
            % reject the same epochs from other erpType epochs as well
            % excluding the EOG and ECG
            rejectFields = {'RAW'; 'EOG'; 'ECG'};
            stimTypes = {'target'; 'distr'; 'std'};
            listOfErpTypes = fieldnames(epochs.(filteringType));
            [erpTypes, indicesOfRemaining] = setdiff(listOfErpTypes, rejectFields);
            for erp = 1 : length(erpTypes)
                for stim = 1 : length(stimTypes)
                    epochs.(filteringType).(erpTypes{erp}).(stimTypes{stim}) = ...
                        pre_rejectEpochsBasedOnFASTER(epochs.(filteringType).(erpTypes{erp}).(stimTypes{stim}), artifactIndices_FASTER.(stimTypes{stim}), stimTypes{stim}, erpTypes{erp}, handles.parameters, handles);
                end
            end
                                                  
        else
            
            disp('     skipping FASTER Artifact Rejection, are you sure?')
            artifactIndices_FASTER = [];
            artifactIndices_FASTER_distr = [];
            artifactIndices_FASTER_std = [];            
            
        end
        close all
        disp(' ')
      
    %% DENOISING (SINGLE-TRIAL ERPs)
    
        filteringType = 'bandpass';
        erpTypes = {'EP_DEN'};
        stimTypes = {'target'; 'distr'; 'std'};         
        
        for erpType = 1 : length(erpTypes)
            
            %% PRE-PROCESS the data for the DENOISING. 
            % i.e. the EP DENOISING needs the epochs to be concatenated together
            for stim = 1 : length(stimTypes)
                epochs_concan.(filteringType).(erpTypes{erpType}).(stimTypes{stim}) = ...
                    pre_concatenateEpochs(epochs.(filteringType).(erpTypes{erpType}).(stimTypes{stim}));
            end          
            
            %% EP_den_auto 
            % ---------------------------------------------------------------------------------------
            % Ahmadi M, Quian Quiroga R. 2013. 
            % Automatic denoising of single-trial evoked potentials. 
            % NeuroImage 66:672–680. http://dx.doi.org/10.1016/j.neuroimage.2012.10.062.
            % Code: http://www2.le.ac.uk/centres/csn/software/ep_den                                       
            filteringTypeOut = 'EP';
                
            for stim = 1 : length(stimTypes)                
                
                % There is only one parameter for the EP DENOISING (that we
                % use)  that is the scale parameter, the larger that is, the
                % smoother is the epoch, you might wanna modify this
                % parameter for different ERP Types, CNV might use a
                % larger value than post-onset components for example
                % additionally you can switch from 'NZT' to 'Neigh' method
                scale = handles.parameters.ep_den.scales_postStim;
                
                if handles.parameters.ep_den.denoiseWithEP == 1
                    if stim == 1
                        disp('    EP Denoising the ERP epochs (P300 component)')
                    end
                    disp(['         .. ', stimTypes{stim}])
                    epochs_concan.(filteringTypeOut).(erpTypes{erpType}).(stimTypes{stim}) = ...
                        denoise_ep_den_auto_Wrapper(epochs_concan.(filteringType).(erpTypes{erpType}).(stimTypes{stim}), ...
                        scale, handles.parameters, handles);

                    % Additional post-smoothing (Jongsma et al. used 3-point moving
                    % average for example)
                else
                    if stim == 1
                        disp('    skipping EP Denoising for the ERP epochs (P300 component)')                                
                    end
                    epochs_concan.(filteringType).(erpTypes{erpType}).(stimTypes{stim}) = [];
                end

            end               
                
            
            %% STEP / N1 measure, 
            % ----------------------------------------------------------------------------------------
            % Hu L, Mouraux A, Hu Y, Iannetti GD. 2010. 
            % A novel approach for enhancing the signal-to-noise ratio and detecting automatically event-related potentials (ERPs) in single trials. 
            % NeuroImage 50:99–111. http://dx.doi.org/10.1016/j.neuroimage.2009.12.010.
            % Code: http://iannettilab.webnode.com/n1measure/
            handles.parameters.ep_den.denoiseWithSTEP = 0;
            filteringTypeOut = 'STEP';
            for stim = 1 : length(stimTypes)               
            
                if handles.parameters.ep_den.denoiseWithSTEP == 1
                
                    if stim == 1
                        disp('    STEP Denoising the ERP epochs (N1 Component)')
                    end
                    disp(['         .. ', stimTypes{stim}])
                    % implement later if needed
                    %{
                    epochs_concan.(filteringType).(erpTypes{erpType}).(stimTypes{stim}) = ...
                        denoise_n1measure_wrapper(epochs_concan.(filteringType).(erpTypes{erpType}).(stimTypes{stim}), ...
                        handles.parameters.ep_den.scales_preStim, handles.parameters, handles);                    
                    %}
                else
                    if stim == 1
                        disp('    skipping STEP Denoising the ERP epochs (N1 Component)')                                
                    end
                    % epochs_concan.(filteringType).(erpTypes{erpType}).(stimTypes{stim}) = [];
                end
            end         

            %% Deconcatenate back to cells of individual epochs
            rejectFields = {'bandpass'};            
            listOfFilterTypes = fieldnames(epochs_concan); % might need to change if not all the denoising method require this concatenation
            [filterTypes, indicesOfRemaining] = setdiff(listOfFilterTypes, rejectFields);
            
            for filt = 1 : length(filterTypes)
                for stim = 1 : length(stimTypes)                                
                    epochs.(filterTypes{filt}).(erpTypes{erpType}).(stimTypes{stim}) ...
                        = pre_deconcatenateEpochs(epochs_concan.(filterTypes{filt}).(erpTypes{erpType}).(stimTypes{stim}), ...
                        handles.parameters, handles);   
                end
            end
            
        end % end of erpTypes (ERP, CNV, ALPHA, ETC.)          
                

    %% Analyze the denoised/filtered epochs, i.e. get the ERP component amplitudes and latencies            

        % Get amplitudes, latencies, etc.
        disp(' '); disp('    Analyzing the ERP components')
        
        filterTypes = fieldnames(epochs);        
        stimTypes = fieldnames(epochs.(filterTypes{1}).(erpTypes{1}));
        
        for filt = 1 : length(filterTypes)
            disp(['     . ', filterTypes{filt}])
            erpTypes = fieldnames(epochs.(filterTypes{filt}));
            
            for erp = 1 : length(erpTypes)
                disp(['      .. ', erpTypes{erp}])
                
                for stim = 1 : length(stimTypes)
                    
                    if stim == 1; fprintf('       ... '); end
                    fprintf([' / ', stimTypes{stim}])
                
                    % Use subfunction
                    ERP_components.(filterTypes{filt}).(erpTypes{erpType}).(stimTypes{stim}) = ...
                        analyze_getERPcomponents(epochs.(filterTypes{filt}).(erpTypes{erpType}).(stimTypes{stim}), ...
                            filterTypes{filt}, handles.parameters.oddballTask.timeWindows, handles.parameters, handles);
                   
                    if stim == length(stimTypes); fprintf(['\n']); end
                        
                end % end of stimulus types (target, distracter, standard)
                
            end % end of ERP types (ERP, CNV, etc.)
            
        end % end of filtering types (Bandpass, EP, etc.)       

    %% TIME-FREQUENCY ANALYSIS FOR THE EPOCHS
    disp(' ')
    if handles.parameters.timeFreq.computeTF == 1
        
        % ~ 805 sec with per file
        disp('    Time-Frequency Analysis (Morlet, CWT)')          

            filteringType = 'bandpass';
            erpType = 'GENERAL';
            [timeFreqEpochs, timeFreq.target, timeFreq.distr, timeFreq.std, TF_allEpochs, TF_derivedMeasures] = ...
                    pre_waveletWrapperForAllStimuli(epochs.(filteringType).(erpType).target, artifactIndices_FASTER.target, ... % Target
                                                    epochs.(filteringType).(erpType).distr, artifactIndices_FASTER.distr, ... % Distracter
                                                    epochs.(filteringType).(erpType).std, artifactIndices_FASTER.std, ... % Standard
                                                    alpha.amplit_gravity, handles.parameters, handles);
    else
        
        % ~ 240 sec without per file
        disp('    Skipping Time-Frequency Analysis (Morlet, CWT)')  
        % assign something that the DISK SAVING function won't crash
        timeFreqEpochs = [];
        timeFreq = [];
        TF_allEpochs = [];
        TF_derivedMeasures = [];
        
    end
    
    %% PHASIC CARDIAC ANALYSIS (PCR) for the EPOCHS
    disp(' '); disp('    Phasic Cardiac Analysis (PCR, Kardia)  (placeholder)')
    
        % Use KARDIA
        % http://sourceforge.net/projects/mykardia/
        % http://dx.doi.org/10.1016/j.cmpb.2009.10.002
        filteringType = 'bandpass';
        erpType = 'ECG';
        stimTypes = {'target'; 'distr'; 'std'};
        for stim = 1 : length(stimTypes)
            heart.PCR.(stimTypes{stim}) = analyze_phasicCardiacResponse(epochs.(filteringType).(erpType).(stimTypes{stim}), handles.parameters, handles);
        end 
        
    %% EVENT-RELATED EOG
    disp(' '); disp('    Event-related EOG (placeholder)')
    
        % not probably anything relevant especially as only one EOG lead is
        % used, but you have the epochs here if you want to further analyze
        % them
        filteringType = 'bandpass';
        erpType = 'EOG';
        stimTypes = {'target'; 'distr'; 'std'};
        for stim = 1 : length(stimTypes)
            EOG.event.(stimTypes{stim}) = analyze_eventRelatedEOG(epochs.(filteringType).(erpType).(stimTypes{stim}), handles.parameters, handles);
        end  

    %% Package to structures and SAVE TO DISK
    disp(' '); disp('    SAVE computed variables to disk')
        
        % you can discard some data, and make sure that you only the needed
        % data is saved to disk, or if you just want to optimize the use of
        % disk space or something    
        process_assignOutputsFromPROCESS(epochs, ERP_components, alpha, powers, SEM, heart, EOG, fractalAnalysis, timeFreqEpochs, timeFreq, TF_allEpochs, TF_derivedMeasures, handles.parameters, handles)    
      
                          
    disp('+++++         Processing of the file complete')
    disp(' ')    
    close all