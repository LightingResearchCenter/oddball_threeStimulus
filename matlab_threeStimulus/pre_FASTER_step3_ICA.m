function [EEG, indelec_st3, zs_st3, num_pca, activData, blinkData] = pre_FASTER_step3_ICA(EEGmatrix, EEG, k_value, ica_chans, chans_to_interp, lpf_band, blinkCh, epochLength, parameters)
            
    % quick'n'dirty
    EEG.dataIN = EEG.data;
    EEG.data = EEGmatrix';

    num_pca = min(floor(sqrt(size(EEG.data(:,:),2) / k_value)),(size(EEG.data,1) - length(chans_to_interp) - 1));
    num_pca = min(num_pca,length(setdiff(ica_chans,chans_to_interp)));

    if num_pca == 0
       warning(['num_pca = ', num2str(num_pca)]) 
       num_pca = 1;
    end

    try

        % Original infomax implementation (slow)                
        [EEG.icaweights, EEG.icasphere, compvars, bias, signs, lrates, EEG.icaact] = runica(EEGmatrix', 'extended', 1, 'pca', num_pca, 'verbose', 'off');
        unmixing_matrix = EEG.icaweights*EEG.icasphere;

        % We could use FastICA instead, suggested also in the discussion
        % of FASTER, http://research.ics.aalto.fi/ica/fastica/
        %{
        [A, unmixing_matrix] = fastica(EEGmatrix', 'lastEig', num_pca, 'verbose', 'off', 'displayMode', 'off'); % gives only the estimated mixing matrix A and the separating matrix W.

            %size(EEG.icaweights), % number of PCAs x number of channels
            %size(EEG.icasphere), % number of PCAs x number of ch
            %size(unmixing_matrix)

            EEG.icaweights = A'; % is this correct?
            %EEG.icasphere ?

            % how to define the number of PCAs, and EXTENDED?
            % check that 'lastEig' is the same as above for runica


        % compute ICA activation waveforms = weights*sphere*(data-meandata)
        % Usage: >> [activations] = icaact(data,weights,datamean);
        %}

        % EEGLAB variables, see e.g. http://sccn.ucsd.edu/wiki/A05:_Data_Structures                
            EEG.icachansind = ica_chans;
            EEG.trials = length(EEGmatrix) / epochLength;
            EEG.pnts = epochLength;                

        EEG.srate = parameters.EEG.srate;
        EEG.icaact = icaact(EEGmatrix', unmixing_matrix);
        EEG.icawinv = pinv(EEG.icaweights*EEG.icasphere); % http://sccn.ucsd.edu/pipermail/eeglablist/2009/002907.html                

        % size(EEG.icaact) % number of PCAs x dataSamples


    catch err
        err
        error('improve error catching!!')
    end

    % after the ICA routine we have the EEG data as 2-dimensional
    % matrix, and we need to to separate the epochs to the third
    % dimension
    [list_properties, activData, blinkData] = component_properties_mod(EEG, blinkCh,lpf_band);
    rejection_options.measure=ones(1,size(list_properties,2)); % values of the statistical parameters (see flow chart)
    rejection_options.z = parameters.artifacts.FASTER_zThreshold * ones(1,size(list_properties,2)); % Z-score threshold
    [indelec_st3, zs_st3] = min_z_mod(list_properties,rejection_options); % rejected components

    step3_linearIndices = find(indelec_st3);
    if ~isempty(step3_linearIndices)
        disp(['          - subtracting ICA artifacts, found ', num2str(length(step3_linearIndices)), ' artifacted ICA activation channels']) 
        for i = 1 : length(step3_linearIndices)
            for ch = 1 : ica_chans
                EEG.dataIN(:,ch) = EEG.icaact(:, step3_linearIndices(i));
            end
        end
    else
        disp(['          - No ICA artifacts found'])             
    end

    % quick'n'dirty
    EEG.data = EEG.dataIN;
    