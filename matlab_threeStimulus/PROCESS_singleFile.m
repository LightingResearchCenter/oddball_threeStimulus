function [epochs_target, epochs_distracter, epochs_standard, ...
          epochsEP_target, epochsEP_distracter, epochsEP_standard, ...
          analyzed_target, analyzed_distracter, analyzed_standard, ...
          analyzed_aux, analyzed_extraSensors, analyzed_fractal, analyzed_TF, ...
          dataMatrix_filtGeneral, alpha, powers, handles] = ...
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

    %% PRE-PROCESS THE DATA            
    [dataMatrix_filtGeneral, firstBandpassFilteredMatrix, artifactNaN_indices, dataMatrix_filtAlpha, dataMatrix_filt, dataMatrix_filt_CNV] = process_preProcessFiltering(dataMatrixIn, handles);
    
            
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
        [alpha, powers, amplitSpectrum, PSD, SEM, heart, fractalAnalysis] =  process_timeSeriesEEG(dataMatrix_filtGeneral(:,EEGind(1):EEGind(2)), ... % EEG
                                                                            dataMatrix_filtGeneral(:,EOGind(1):EOGind(2)), ... % EOG 
                                                                            dataMatrixIn(:,EOGind_raw(1):EOGind_raw(2)), ... % EOG RAW
                                                                            dataMatrix_filtGeneral(:,ECGind(1):ECGind(2)), ... % ECG
                                                                            dataMatrixIn(:,ECGind_raw(1):ECGind_raw(2)), ... % ECG RAW
                                                                            triggers, handles.style, handles.parameters, handles);                                                                           
        close all
                       
    %% Epoch the EEG
    
        % i.e. Split into Oddball, Distracter and Standard Event-Related Potentials (ERPs)            
        disp('    Split recording to ERP epochs')              
        
        % process_epochingToERPs(dataMatrix_filt, dataMatrix_filtGeneral, dataMatrix_filt_CNV, dataMatrix_filtAlpha, triggers, alpha, handles)
            
            disp('     ERP FILTERED')  
            [epochs_filt, epochs_distr_filt, epochs_std_filt, epoch_indices] = pre_epochToERPs(dataMatrix_filt(:,:), triggers, [], alpha, handles.parameters, 'filt', handles);    

            disp('     RAW INPUT')  
            [epochs_rawInput, epochs_distr_rawInput, epochs_std_rawInput, ~] = pre_epochToERPs(dataMatrixIn(:,:), triggers, epoch_indices, alpha, handles.parameters, 'rawInput', handles);       
                        
            disp('     GENERAL')  
            [epochs_raw, epochs_distr_raw, epochs_std_raw, ~] = pre_epochToERPs(dataMatrix_filtGeneral(:,:), triggers, epoch_indices, alpha, handles.parameters, 'raw', handles);
            
            disp('     CNV FILTERED')  
            [epochs_CNV_filt, epochs_CNV_distr_filt, epochs_CNV_std_filt, ~] = pre_epochToERPs(dataMatrix_filt_CNV(:,:), triggers, epoch_indices, alpha, handles.parameters, 'CNV', handles);
            
            disp('     ALPHA FILTERED')  
            [epochs_filtAlpha, epochs_distr_filtAlpha, epochs_std_filtAlpha, ~] = pre_epochToERPs(dataMatrix_filtAlpha(:,:), triggers, epoch_indices, alpha, handles.parameters, 'Alpha', handles);
            
            disp('     ECG')
            [epochs_ECG, epochs_distr_ECG, epochs_std_ECG, ~] = pre_epochToERPs(dataMatrixIn(:,handles.parameters.EEG.nrOfChannels+2:handles.parameters.EEG.nrOfChannels+2), triggers, epoch_indices, alpha, handles.parameters, 'ECG', handles);
            
            disp('     EOG')
            [epochs_EOG, epochs_distr_EOG, epochs_std_EOG, ~] = pre_epochToERPs(dataMatrix_filtGeneral(:,handles.parameters.EEG.nrOfChannels+1:handles.parameters.EEG.nrOfChannels+1), triggers, epoch_indices, alpha, handles.parameters, 'ECG', handles);
            
                % maybe change the ALPHA to actual time-frequency
                % presentation, or just add theta also :
                
                % SEE Polich (2007), pg. 240-241
                
                % Polich J. 2007. 
                % Updating P300: An integrative theory of P3a and P3b. 
                % Clinical Neurophysiology 118:2128–2148. 
                % http://dx.doi.org/10.1016/j.clinph.2007.04.019.            

            % save EOG and ECG
            %{
            EOG = dataMatrix_filt(:,7-offset);
            EOG_CNV = dataMatrix_filt_CNV(:,7-offset);
            ECG = dataMatrix_filt(:,8-offset);
            ECG_CNV = dataMatrix_filt_CNV(:,8-offset);
            %}
        
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
            disp('    FASTER Artifact Rejection')           
            [epochs_raw, epochs_distr_raw, epochs_std_raw, artifactIndices_FASTER, artifactIndices_FASTER_distr, artifactIndices_FASTER_std]...
                = pre_FASTER_forAllStimuli(epochs_raw, epochs_distr_raw, epochs_std_raw, ...
                                           epochs_rawInput, epochs_distr_rawInput, epochs_std_rawInput, ...
                                           epochs_EOG, epochs_distr_EOG, epochs_std_EOG, ...
                                           handles.parameters, handles);
        else
            
            disp('     skipping FASTER Artifact Rejection, are you sure?')
            artifactIndices_FASTER = [];
            artifactIndices_FASTER_distr = [];
            artifactIndices_FASTER_std = [];            
            
        end
        close all
        

    %% TIME-FREQUENCY ANALYSIS FOR THE EPOCHS    
    disp('    Time-Frequency Analysis (Morlet, CWT)')  

        [timeFreqEpochs, timeFreq_target, timeFreq_distr, timeFreq_std] = ...
                pre_waveletWrapperForAllStimuli(epochs_raw, artifactIndices_FASTER, epochs_distr_raw, artifactIndices_FASTER_distr, epochs_std_raw, artifactIndices_FASTER_std, ...
                alpha.amplit_gravity, handles.parameters, handles);
    
    %% PHASIC CARDIAC ANALYSIS (PCR) for the EPOCHS
    disp('    Phasic Cardiac Analysis (PCR, Kardia)  (placeholder)')
    
        % Use KARDIA
        % http://sourceforge.net/projects/mykardia/
        % http://dx.doi.org/10.1016/j.cmpb.2009.10.002
        heart.PCR.target = analyze_phasicCardiacResponse(epochs_ECG, handles.parameters, handles);
        heart.PCR.distracter = analyze_phasicCardiacResponse(epochs_distr_ECG, handles.parameters, handles);
        heart.PCR.standard = analyze_phasicCardiacResponse(epochs_std_ECG, handles.parameters, handles);
        
    %% EVENT-RELATED EOG
    disp('    Event-related EOG (placeholder)')
    
        % not probably anything relevant especially as only one EOG lead is
        % used, but you have the epochs here if you want to further analyze
        % them
        EOG.event.target = analyze_eventRelatedEOG(epochs_EOG, handles.parameters, handles);
        EOG.event.distracter = analyze_eventRelatedEOG(epochs_distr_EOG, handles.parameters, handles);
        EOG.event.standard = analyze_eventRelatedEOG(epochs_std_EOG, handles.parameters, handles);
        
            
    %% PRE-PROCESS the data for the DENOISING
        disp(' ')    
                       
        % EP_den auto requires all the epochs to be concatenated into a
        % single vector (single vector per channel)
        epochs_concan = pre_concatenateEpochs(epochs_filt, handles.parameters, handles);        
        epochs_concan_CNV = pre_concatenateEpochs(epochs_CNV_filt, handles.parameters, handles);
        
        epochs_concan_distr = pre_concatenateEpochs(epochs_distr_filt, handles.parameters, handles);        
        epochs_concan_CNV_distr = pre_concatenateEpochs(epochs_CNV_distr_filt, handles.parameters, handles);
        
        epochs_concan_std = pre_concatenateEpochs(epochs_std_filt, handles.parameters, handles);        
        epochs_concan_CNV_std = pre_concatenateEpochs(epochs_CNV_std_filt, handles.parameters, handles);        
       

    %% DENOISE the EPOCHS to obtain single-trial ERPs without too much noise               

        % EP_den_auto 
        % ---------------------------------------------------------------------------------------
        % Ahmadi M, Quian Quiroga R. 2013. 
        % Automatic denoising of single-trial evoked potentials. 
        % NeuroImage 66:672–680. http://dx.doi.org/10.1016/j.neuroimage.2012.10.062.
        % Code: http://www2.le.ac.uk/centres/csn/software/ep_den
        handles.parameters.ep_den.denoiseWithEP = 1;
        if handles.parameters.ep_den.denoiseWithEP == 1
        
            disp('    EP Denoising the ERP epochs')

            % ODDBALL / TARGET
            
                disp('         .. ODDBALL / TARGET')
                % 1st to obtain the CNV            
                epochs_ep_CNV = denoise_ep_den_auto_Wrapper(epochs_concan_CNV, handles.parameters.ep_den.scales_preStim, handles.parameters, handles);

                % 2nd to obtain N2,P3 
                epochs_ep = denoise_ep_den_auto_Wrapper(epochs_concan, handles.parameters.ep_den.scales_postStim, handles.parameters, handles);
                
            % DISTRACTER / non-TARGET
            
                disp('         .. DISTRACTER / NON-TARGET')
                % 1st to obtain the CNV            
                epochs_ep_CNV_distr = denoise_ep_den_auto_Wrapper(epochs_concan_CNV_distr, handles.parameters.ep_den.scales_preStim, handles.parameters, handles);

                % 2nd to obtain N2,P3 
                epochs_ep_distr = denoise_ep_den_auto_Wrapper(epochs_concan_distr, handles.parameters.ep_den.scales_postStim, handles.parameters, handles);
                
            % STANDARD
            
                disp('         .. STANDARD')
                % 1st to obtain the CNV            
                epochs_ep_CNV_std = denoise_ep_den_auto_Wrapper(epochs_concan_CNV_std, handles.parameters.ep_den.scales_preStim, handles.parameters, handles);

                % 2nd to obtain N2,P3 
                epochs_ep_std = denoise_ep_den_auto_Wrapper(epochs_concan_std, handles.parameters.ep_den.scales_postStim, handles.parameters, handles);
        else
            
            disp('    skipping EP Denoising the ERP epochs')
            epochs_ep_CNV = [];
            epochs_ep = [];
            epochs_ep_CNV_distr = [];
            epochs_ep_distr = [];
            epochs_ep_CNV_std = [];
            epochs_ep_std = [];
            
        end
            

        % STEP / N1 measure, 
        % ----------------------------------------------------------------------------------------
        % Hu L, Mouraux A, Hu Y, Iannetti GD. 2010. 
        % A novel approach for enhancing the signal-to-noise ratio and detecting automatically event-related potentials (ERPs) in single trials. 
        % NeuroImage 50:99–111. http://dx.doi.org/10.1016/j.neuroimage.2009.12.010.
        % Code: http://iannettilab.webnode.com/n1measure/
        % epochs_N1 = denoise_n1measure_wrapper(epochs, handles.parameters, handles);

            % Might be hard to distinguish N1 from N2, and in the
            % papers of Jongsma et al., no N1 was analyzed so this
            % analysis is not done at the moment.

        % Additional post-smoothing (Jongsma et al. used 3-point moving
        % average if needed        

            % something if wanted

        % Deconcatenate back to cells of individual epochs
        epochs_ep = pre_deconcatenateEpochs(epochs_ep, handles.parameters, handles);
        epochs_ep_CNV = pre_deconcatenateEpochs(epochs_ep_CNV, handles.parameters, handles);

        epochs_ep_distr = pre_deconcatenateEpochs(epochs_ep_distr, handles.parameters, handles);
        epochs_ep_CNV_distr = pre_deconcatenateEpochs(epochs_ep_CNV_distr, handles.parameters, handles);
        
        epochs_ep_std = pre_deconcatenateEpochs(epochs_ep_std, handles.parameters, handles);
        epochs_ep_CNV_std = pre_deconcatenateEpochs(epochs_ep_CNV_std, handles.parameters, handles);
        
        % Trim the epochs so that the pre-onset baseline is removed,
        % input is now back in Cell (not concatenated)
            % epochs_ep = pre_removeBaseline(epochs_ep, handles.parameters, handles);
            % epochs_raw = pre_removeBaseline(epochs, handles.parameters, handles);

    %% Analyze the denoised/filtered epochs, i.e. get the ERP component amplitudes and latencies            

        % Get amplitudes, latencies, etc.

            % YOU MIGHT WANT TO ADD A SWITCH INSTEAD analyze_ at some point
            % to select whether you wanna do EP or FILTERED.. this is just
            % a quick'n'dirty fix
            
            % i.e. too many inputs going in to the subfunction
        
        % ODDBALL / TARGET
        disp('    Analyzing the ERP components')
        disp('            Oddballs/Targets')
        [ERP_components, epochs_ERP, epochs_ERP_CNV, epochs_ERP_raw, epochs_ERP_filt, epochs_ERP_CNV_filt, handles.parameters.oddballTask.timeVector] = ...
                analyze_getERPcomponents(epochs_filt, epochs_CNV_filt, epochs_raw, epochs_filt, epochs_CNV_filt, 'filt', handles.parameters.oddballTask.timeWindows, handles.parameters, handles);
            
        [ERP_components_EP, ~, ~, ~, ~, ~, ~] = ...
                analyze_getERPcomponents(epochs_ep, epochs_ep_CNV, epochs_raw, epochs_filt, epochs_CNV_filt, 'EP', handles.parameters.oddballTask.timeWindows, handles.parameters, handles);
            
        % DISTRACTER / non-TARGET
        disp('            Distracters/non-Targets')
        [ERP_components_distr, epochs_ERP_distr, epochs_ERP_CNV_distr, epochs_ERP_raw_distr, epochs_ERP_filt_distr, epochs_ERP_CNV_filt_distr, handles.parameters.oddballTask.timeVector] = ...
                analyze_getERPcomponents(epochs_distr_filt, epochs_CNV_distr_filt, epochs_distr_raw, epochs_distr_filt, epochs_CNV_distr_filt, 'filt', handles.parameters.oddballTask.timeWindows, handles.parameters, handles);
            
        [ERP_components_EP_distr, ~, ~, ~, ~, ~, ~] = ...
                analyze_getERPcomponents(epochs_ep_distr, epochs_ep_CNV_distr, epochs_distr_raw, epochs_distr_filt, epochs_CNV_distr_filt, 'EP', handles.parameters.oddballTask.timeWindows, handles.parameters, handles);
        
        % STANDARD
        disp('            Standards')
        [ERP_components_std, epochs_ERP_std, epochs_ERP_CNV_std, epochs_ERP_raw_std, epochs_ERP_filt_std, epochs_ERP_CNV_filt_std, handles.parameters.oddballTask.timeVector] = ...
                analyze_getERPcomponents(epochs_std_filt, epochs_CNV_std_filt, epochs_std_raw, epochs_std_filt, epochs_CNV_std_filt, 'filt', handles.parameters.oddballTask.timeWindows, handles.parameters, handles);
            
        [ERP_components_EP_std, ~, ~, ~, ~, ~, ~] = ...
                analyze_getERPcomponents(epochs_ep_std, epochs_ep_CNV_std, epochs_std_raw, epochs_std_filt, epochs_CNV_std_filt, 'EP', handles.parameters.oddballTask.timeWindows, handles.parameters, handles);
        
                    % ADD ALPHA


    %% Package to structures 
    
        % three main types
        
            % * "classical" bandpass filtered ERPs
            % * "EP denoised"
            % The scalars (amplitudes, latencies) extracted from ERPs
    
    
        %% EPOCHS IN        
            
            epochs_target.ERP = epochs_filt;
            epochs_target.CNV = epochs_CNV_filt;
            epochs_target.General = epochs_raw;
            epochs_target.Alpha = epochs_filtAlpha;
            
            epochs_distracter.ERP = epochs_distr_filt;
            epochs_distracter.CNV = epochs_CNV_distr_filt;
            epochs_distracter.General = epochs_distr_raw;
            epochs_distracter.Alpha = epochs_distr_filtAlpha;            
            
            epochs_standard.ERP = epochs_std_filt;
            epochs_standard.CNV = epochs_CNV_std_filt;
            epochs_standard.General = epochs_std_raw;
            epochs_standard.Alpha = epochs_std_filtAlpha;           
            
            
        
        %% EP denoised EPOCHS        
       
            epochsEP_target.ERP = epochs_ep;
            epochsEP_target.CNV = epochs_ep_CNV;
            % general
            % alpha
            
            epochsEP_distracter.ERP = epochs_ep_distr;
            epochsEP_distracter.CNV = epochs_ep_CNV_distr;            
            % general
            % alpha
            
            epochsEP_standard.ERP = epochs_ep_std;
            epochsEP_standard.CNV = epochs_ep_CNV_std;
            % general 
            % alpha
    
        %% Analyzed 
        
            % target/oddball
            analyzed_target.ERPcomponents.filt = ERP_components;
            analyzed_target.ERPcomponents.EP = ERP_components_EP;
            analyzed_target.ERP = epochs_ERP;
            analyzed_target.CNV = epochs_ERP_CNV;
            analyzed_target.ERP_general = epochs_ERP_raw;
            analyzed_target.ERP_filt = epochs_ERP_filt;
            analyzed_target.CNV_filt = epochs_ERP_CNV_filt;
        
            % distracter / non-target
            analyzed_distracter.ERPcomponents.filt = ERP_components_distr;
            analyzed_distracter.ERPcomponents.EP = ERP_components_EP_distr;
            analyzed_distracter.ERP = epochs_ERP_distr;
            analyzed_distracter.CNV = epochs_ERP_CNV_distr;
            analyzed_distracter.ERP_general = epochs_ERP_raw_distr;
            analyzed_distracter.ERP_filt = epochs_ERP_filt_distr;
            analyzed_distracter.CNV_filt = epochs_ERP_CNV_filt_distr;
        
            % standard
            analyzed_standard.ERPcomponents.filt = ERP_components_std;
            analyzed_standard.ERPcomponents.EP = ERP_components_EP_std;
            analyzed_standard.ERP = epochs_ERP_std;
            analyzed_standard.CNV = epochs_ERP_CNV_std;
            analyzed_standard.ERP_general = epochs_ERP_raw_std;
            analyzed_standard.ERP_filt = epochs_ERP_filt_std;
            analyzed_standard.CNV_filt = epochs_ERP_CNV_filt_std;
        
            % auxiliary measures (scalars)
            analyzed_aux.IAF_amplitGravity = alpha.amplit_gravity;
            
            for i = 1 : length(handles.parameters.powerAnalysis.eegBins.freqs)
                fieldName = powers{i}.label;
                analyzed_aux.Amplit.(fieldName) = powers{i}.ps.powerData;                
                analyzed_aux.PSD.(fieldName) = powers{i}.PSD.powerData;
            end
            
            %analyzed_aux
            %analyzed_aux.Amplit
            
            % EOG and ECG Analysis
            analyzed_extraSensors.SEM = SEM;
            analyzed_extraSensors.heart = heart;
               
            % Fractal analysis for individual channels
            analyzed_fractal = fractalAnalysis;
            
            % Wavelet time-frequency analysis for the epochs, 
           % i.e. to find ERD
            analyzed_TF.target = timeFreq_target;
            analyzed_TF.distr = timeFreq_distr;
            analyzed_TF.std = timeFreq_std;
            
                          
    disp('+++++         Processing of the file complete')
    disp(' ')
    
    close all