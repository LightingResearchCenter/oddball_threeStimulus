function [rejectedEpochs_target, rejectedEpochs_distr, rejectedEpochs_std, artifactIndices_FASTER_target, artifactIndices_FASTER_distr, artifactIndices_FASTER_std] = ...
    pre_FASTER_forAllStimuli(rejectOn_target, rejectOn_distr, rejectOn_std, ... % reject
    reference_target, reference_distr, reference_std, ... % reference
    EOG_target, EOG_distr, EOG_std, ... % EOG
    ECG_target, ECG_distr, ECG_std, ... % ECG
    filteringType, rejectBasedOn, referenceType, ... % strings 
    debugOn, ... % boolean flags
    parameters, handles) % parameters
    
    localCall = 0;            

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempFASTER_allStimuli.mat';
        if nargin == 0
            load('debugPath.mat')
            load(fullfile(path.debugMATs, debugMatFileName))
            close all
            localCall = 1;
        else
            if handles.flags.saveDebugMATs == 1
                path = handles.path;
                save('debugPath.mat', 'path')
                save(fullfile(path.debugMATs, debugMatFileName))
            end
        end 
    end
    
    %{
    handles.parameters.artifacts.FASTER_zThreshold_step2 = 1.49;
    handles.parameters.artifacts.FASTER_zThreshold_step3 = 1.49;
    handles.parameters.artifacts.FASTER_zThreshold_step4 = 1.49;
    %}

    debugOnlyTarget = 0; % faster development / debugging
    
    disp('      .. Find artifacts (Fixed EEG/EOG, CRAP: Moving Average & Step)')  
        
        % TARGET
        fprintf('         .. TARGET / ')
        
        %% NOTE! (19 Dec 2013)
        % Now we haven't done anything yet for the trials that do not have
        % response (i.e. the TARGET tones with NaN RT (reaction time), see e.g.
        % https://www.researchgate.net/post/Do_we_remove_trials_with_erroneous_RT_answers_in_ERP_analysis
        
        [NaN_indices_EEG_target, NaN_indices_EOG_target, NaN_indices_moving_target, NaN_indices_movingEOG_target, NaN_indices_step_target, ...
            vDiffOutMovWindow_target, vDiffOutMovWindowEOG_target, vDiffOutStep_target, fixedIndices_target] = ...
            pre_artifactFASTER_fixedThresholds_ERPLAB(rejectOn_target, EOG_target, ECG_target, debugOn, 'target', rejectOn_target.samplesPerEpoch, handles.parameters, handles);    
        
        if debugOnlyTarget ~= 1 || localCall == 0
            
            % DISTRACTER
            fprintf('DISTRACTER / ')
            [NaN_indices_EEG_distr, NaN_indices_EOG_distr, NaN_indices_moving_distr, NaN_indices_movingEOG_distr, NaN_indices_step_distr, ...
                vDiffOutMovWindow_distr, vDiffOutMovWindowEOG_distr, vDiffOutStep_distr, fixedIndices_distr] = ...
                pre_artifactFASTER_fixedThresholds_ERPLAB(rejectOn_distr, EOG_distr, ECG_distr, debugOn, 'distracter', rejectOn_distr.samplesPerEpoch, handles.parameters, handles);

            % STANDARD
            fprintf('STANDARD'); fprintf('\n')
            [NaN_indices_EEG_std, NaN_indices_EOG_std, NaN_indices_moving_std, NaN_indices_movingEOG_std, NaN_indices_step_std, ...
                vDiffOutMovWindow_std, vDiffOutMovWindowEOG_std, vDiffOutStep_std, fixedIndices_std] = ...
                pre_artifactFASTER_fixedThresholds_ERPLAB(rejectOn_std, EOG_std, ECG_std, debugOn, 'standard', rejectOn_std.samplesPerEpoch, handles.parameters, handles);
        end

    % concatenate
    epochs_concan_target = pre_concatenateEpochs(rejectOn_target, handles.parameters, handles);
    if debugOnlyTarget ~= 1 || localCall == 0
        epochs_concan_distr = pre_concatenateEpochs(rejectOn_distr, handles.parameters, handles);
        epochs_concan_std = pre_concatenateEpochs(rejectOn_std, handles.parameters, handles);
    end

    % FASTER (3rd party) used with a wrapper funtion
    disp('      .. FASTER algorithm')    
    [epochs_concan_FASTER_target, artifactIndices_FASTER_target] = pre_artifactFASTER_wrapper(epochs_concan_target, fixedIndices_target,...
        NaN_indices_EEG_target, NaN_indices_EOG_target, NaN_indices_moving_target, NaN_indices_movingEOG_target, NaN_indices_step_target, ...
        reference_target, vDiffOutMovWindow_target, vDiffOutMovWindowEOG_target, vDiffOutStep_target, ...
        handles.parameters, 'target', handles);
    
    if debugOnlyTarget ~= 1 || localCall == 0
        [epochs_concan_FASTER_distr, artifactIndices_FASTER_distr] = pre_artifactFASTER_wrapper(epochs_concan_distr, fixedIndices_distr,...
            NaN_indices_EEG_distr, NaN_indices_EOG_distr, NaN_indices_moving_distr, NaN_indices_movingEOG_distr, NaN_indices_step_distr, ...
            reference_distr, vDiffOutMovWindow_distr, vDiffOutMovWindowEOG_distr, vDiffOutStep_distr, ...
            handles.parameters, 'distracter', handles);            

        [epochs_concan_FASTER_std, artifactIndices_FASTER_std] = pre_artifactFASTER_wrapper(epochs_concan_std, fixedIndices_std,...
            NaN_indices_EEG_std, NaN_indices_EOG_std, NaN_indices_moving_std, NaN_indices_movingEOG_std, NaN_indices_step_std, ...
            reference_std, vDiffOutMovWindow_std, vDiffOutMovWindowEOG_std, vDiffOutStep_std, ...
            handles.parameters, 'standard', handles);
    end

    % deconcatenate back
    epochs_deconcan_FASTER = pre_deconcatenateEpochs(epochs_concan_FASTER_target, handles.parameters, handles);
    if debugOnlyTarget ~= 1 || localCall == 0
        epochs_deconcan_FASTER_distr = pre_deconcatenateEpochs(epochs_concan_FASTER_distr, handles.parameters, handles);
        epochs_deconcan_FASTER_std = pre_deconcatenateEpochs(epochs_concan_FASTER_std, handles.parameters, handles);
    end

    % convert the artifacted epochs to NaN vectors, and assign to 
    % different variable names so that the denoising goes
    % okay and you can debug the FASTER step easily if you need
    %{
    rejectOn_target = pre_convertEpochs_to_NanVectors(epochs_deconcan_FASTER, artifactIndices_FASTER);
    rejectOn_distr = pre_convertEpochs_to_NanVectors(epochs_deconcan_FASTER_distr, artifactIndices_FASTER_distr);
    rejectOn_std = pre_convertEpochs_to_NanVectors(epochs_deconcan_FASTER_std, artifactIndices_FASTER_std);
    %}

    % rename to match to upcoming steps
    rejectedEpochs_target = epochs_deconcan_FASTER;
    if debugOnlyTarget ~= 1 || localCall == 0
        rejectedEpochs_distr = epochs_deconcan_FASTER_distr;
        rejectedEpochs_std = epochs_deconcan_FASTER_std;
    end
