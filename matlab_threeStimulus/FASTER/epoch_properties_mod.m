function list_properties = epoch_properties_mod(EEG,EEG_joined,eeg_chans,inputType,epochLength)
    
    list_properties = [];
    measure = 1;

    % Added the input type switch so you can push normal MATLAB matrices
    % instead of having the data in EEGLab format, Petteri Teikari
    if strcmp(inputType, 'EEGLAB')
        if length(size(EEG.data)) < 3
            fprintf('Not epoched.\n');
            return;
        end
    elseif strcmp(inputType, 'matrix')
          
    else
        inputType
        error('Problem with the inputType, a typo!')
    end
    
    % For this function, there should be no NaNs coming in    
    noOfNaNs = sum(sum(sum(isnan(EEG.data))));
    if noOfNaNs > 0
        warning(['There should not be any NaNs (n = ', num2str(noOfNaNs), ') coming to this FASTER function, check what is wrong'])
    end

    % get channel means, for 4 electrodes, you should get 4 values
    EEGmat = EEG_joined(eeg_chans,:);
    try
        means = mean(EEG_joined(eeg_chans,:),2);
    catch err
       err
       which -all mean
       error('Maybe a conflict with BioSig''s mean function? Remove extra from path')
    end

    % Use a bit non-elegant way to duplicate all the steps for EEGLAB and
    % 'matrix' input type

    %% 1 Epoch's mean deviation from channel means.
    for u = 1:size(EEG.data,3)
        list_properties(u,measure) = mean(abs(squeeze(mean(EEG.data(eeg_chans,:,u),2)) - means));
    end    
    measure = measure + 1;

    %% 2 Epoch variance
    list_properties(:,measure) = mean(squeeze(var(EEG.data(eeg_chans,:,:),0,2)));
    measure = measure + 1;

    %% 3 Max amplitude difference
    for t = eeg_chans
        for u = 1:size(EEG.data,3)
            ampdiffs(t,u) = max(EEG.data(t,:,u)) - min(EEG.data(t,:,u));
        end
    end
    list_properties(:,measure) = mean(ampdiffs,1);
    measure = measure + 1;

    for v = 1:size(list_properties,2)
        list_properties(:,v) = list_properties(:,v) - median(list_properties(:,v));
    end



        

        