function [EEG, indelec_st3, zs_st3, num_pca, activData, blinkData] = pre_FASTER_step3_ICA(EEGmatrix, EEG, k_value, ica_chans, chans_to_interp, lpf_band, blinkCh, epochLength, parameters, handles)
            
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempFASTER_ICA.mat';
        if nargin == 0
            load('debugPath.mat')
            load(fullfile(path.debugMATs, debugMatFileName))
            
        else
            if handles.flags.saveDebugMATs == 1                
                % do not save for standard tone as there are so many
                % trials that debugging and developing of this function
                % is so much slower compared to target and distracter
                path = handles.path;
                save('debugPath.mat', 'path')
                save(fullfile(path.debugMATs, debugMatFileName))                            
            end
        end 
    end

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
        
        parameters.artifacts.FASTER_icaMethod = 'runica';
        
        if strcmp(parameters.artifacts.FASTER_icaMethod, 'runica')

            % Original infomax implementation (slow)              
            disp('          ... computing ICA (runica), might take some time (try to switch to fastICA for speed)') 
            [EEG.icaweights, EEG.icasphere, compvars, bias, signs, lrates, EEG.icaact] = runica(EEGmatrix', 'extended', 1, 'pca', num_pca, 'verbose', 'off');
            unmixing_matrix = EEG.icaweights*EEG.icasphere;
            
        elseif strcmp(parameters.artifacts.FASTER_icaMethod, 'fastica')
            
            disp('          ... computing ICA (fastICA), faster than the default runica') 

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
        
        end

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
    
    
    %% EXAMPLE CODE FROM fieldTrip (ft_componentanalysis.m)
    
        % perform the component analysis
        cfg.method = 'fastica'
        if 1 == 2

            switch cfg.method

                case 'fastica'
                % check whether the required low-level toolboxes are installed
                ft_hastoolbox('fastica', 1);       % see http://www.cis.hut.fi/projects/ica/fastica

                % set the number of components to be estimated
                cfg.fastica.numOfIC = cfg.numcomponent;

                try
                  % construct key-value pairs for the optional arguments
                  optarg = ft_cfg2keyval(cfg.fastica);
                  [mixing, unmixing] = fastica(dat, optarg{:});
                catch
                  % the "catch me" syntax is broken on MATLAB74, this fixes it
                  me = lasterror;
                  % give a hopefully instructive error message
                  fprintf(['If you get an out-of-memory in fastica here, and you use fastica 2.5, change fastica.m, line 482: \n' ...
                    'from\n' ...
                    '  if ~isempty(W)                  %% ORIGINAL VERSION\n' ...
                    'to\n' ...
                    '  if ~isempty(W) && nargout ~= 2  %% if nargout == 2, we return [A, W], and NOT ICASIG\n']);
                 % forward original error
                  rethrow(me);
                end

              case 'runica'
                % check whether the required low-level toolboxes are installed
                % see http://www.sccn.ucsd.edu/eeglab
                ft_hastoolbox('eeglab', 1);

                % construct key-value pairs for the optional arguments
                optarg = ft_cfg2keyval(cfg.runica);
                [weights, sphere] = runica(dat, optarg{:});

                % scale the sphering matrix to unit norm
                if strcmp(cfg.normalisesphere, 'yes'),
                  sphere = sphere./norm(sphere);
                end

                unmixing = weights*sphere;
                mixing = [];

            end

            % make sure we have both mixing and unmixing matrices
            % if not, compute (pseudo-)inverse to go from one to the other
            if isempty(unmixing) && ~isempty(mixing)
              if (size(mixing,1)==size(mixing,2))
                unmixing = inv(mixing);
              else
                unmixing = pinv(mixing);
              end
            elseif isempty(mixing) && ~isempty(unmixing)
              if (size(unmixing,1)==size(unmixing,2)) && rank(unmixing)==size(unmixing,1)
                mixing = inv(unmixing);
              else
                mixing = pinv(unmixing);
              end
            elseif isempty(mixing) && isempty(unmixing)
              % this sanity check is needed to catch convergence problems in fastica
              % see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1519
              error('the component unmixing failed');
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % collect the results and construct data structure

            comp = [];
            if isfield(data, 'fsample'), comp.fsample = data.fsample; end
            if isfield(data, 'time'),    comp.time    = data.time;    end

            % make sure we don't return more components than were requested
            % (some methods respect the maxcomponent parameters, others just always
            % return a fixed (i.e., numchans) number of components)
            if size(unmixing,1) > cfg.numcomponent
              unmixing(cfg.numcomponent+1:end,:) = [];
            end

            if size(mixing,2) > cfg.numcomponent
              mixing(:,cfg.numcomponent+1:end) = [];
            end

            % compute the activations in each trial
            for trial=1:Ntrials
              comp.trial{trial} = scale * unmixing * data.trial{trial};
            end

            % store mixing/unmixing matrices in structure
            comp.topo = mixing;
            comp.unmixing = unmixing;

        end

