function EEGcorr = pre_artifactRegressionWrapper(EEG, artifactData, EEGchans, parameters)

    %% EOG/ECG "trend correction"    

    % Requires BioSig
    % http://sourceforge.net/projects/biosig/

    % REGRESS_EOG yields the regression coefficients for 
    % correcting EOG artifacts in EEG recordings.
    % http://biosig.sourceforge.net/help/biosig/t250/regress_eog.html

    % combine EOG/ECG and EEG
    S1 = [EEG artifactData];

    % find the artifacts
    EL = 1:EEGchans;
    OL = EEGchans+1;    

    try
        R = regress_eog(S1, EL, OL);
    catch err
        if strcmp(err.identifier, 'MATLAB:UndefinedFunction')
            error('EOG Rejection not done! Have you installed the BioSig toolbox needed for EOG artifact rejection? http://sourceforge.net/projects/biosig/')
        else
            err
        end                
    end

    % apply the correction, with offset correction
    S2 = [ones(size(S1,1),1),S1] * R.r1;
    % S2b = S1 * R.r0;    % without offset correction
    EEGcorr = S2(:,EL);
 