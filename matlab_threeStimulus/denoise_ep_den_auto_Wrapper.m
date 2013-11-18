function epochs_ep = denoise_ep_den_auto_Wrapper(epochs, scales, parameters, handles)

    debugMatFileName = 'tempDenoiseEP.mat';
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

    % whos
    % parameters
    
    % Ahmadi M, Quian Quiroga R. 2013. 
    % Automatic denoising of single-trial evoked potentials. 
    % NeuroImage 66:672–680. http://dx.doi.org/10.1016/j.neuroimage.2012.10.062.
    % Code: http://www2.le.ac.uk/centres/csn/software/ep_den

    % Summary of the NZT algorithm
    % In summary, the automatic single trial ERPs denoising method proposed here consist of the following steps:  
    % 1. Obtain the average ERP.
    % 2. Do a wavelet decomposition of the average ERP.
    % 3. Do a hard thresholding of the wavelet coefficients using Eq. (6) and the zerotrees procedure.
    % 4. Reconstruct the denoised average ERP using the denoised coefficients.
    % 5. Use the same set of coefficients for the single trial ERPs.

    % assign RTs directly to output
    epochs_ep.RT = epochs.RT;
    
    % also the samples per epoch
    epochs_ep.samplesPerEpoch = epochs.samplesPerEpoch;
    
    % same for indices
    epochs_ep.Indices = epochs.Indices;

    % Loop through the epochs (i.e. the ERPs for oddballs), the channels in
    % other words, as the epochs are concatenated into single vector
    
        den_type = parameters.ep_den.den_type;    
        [rowsIn, colsIn] = size(epochs.ERP);
        
        % EPOCHS
        for i = 1 : parameters.EEG.nrOfChannels
            epochs_ep.ERP(:,i) = epLoopFunction(i, epochs.ERP(:,i), parameters, scales, den_type, epochs.samplesPerEpoch, handles);
        end
        
    
    %% Subfunction for loop
    function epoch = epLoopFunction(i, epochIn, parameters, scales, den_type, samplesPerEpoch, handles)
        
        % Modified the code from EP_den_Auto.m (Batch Folder)    
        samples = samplesPerEpoch;
        sr = parameters.EEG.srate;
        
        % stim defines from where the actual stimulus starts from as the
        % method requires some baseline before the stimulus
        relativeOnset = parameters.oddballTask.ERP_baseline / (parameters.oddballTask.ERP_baseline + parameters.oddballTask.ERP_duration);              
        stim = (sr *  relativeOnset) + 1; % defined in samples
        
        sc = scales;
        plot_type = parameters.ep_den.plot_type;
        
        auto_den_type = parameters.ep_den.auto_den_type;
        
        max_trials = 10;
        max_contour = 32;

        plot_type = 'contour';
        
        %% Denoising
        x = epochIn; % rename 
        sweeps = length(x)/samples;    
        try
            xx=reshape(x,samples,sweeps)';
        catch err
            sweeps
            err
            error('Non-integer nr of sweeps? Problem with the "pre_epochToERPs" definition')
        end                      
            
        av=mean(xx,1);
        
        switch den_type
            
            case 'do_den' 
                
                switch auto_den_type
                    case 'Neigh'
                        try
                            [coeff,denav,den_coeff,y,yo]= Run_Neigh(av,stim,sc);
                        catch err
                            if strcmp(err.identifier, 'MATLAB:TooManyInputs')
                                error('You probably ran the GUI of EP_DEN and now Matlab confuses the two Run_NZT, closing and restarting Matlab at least helps')
                            elseif strcmp(err.identifier, 'MATLAB:nologicalnan')                                
                                error('Why actually? happened with too low high cutoff frequency for bandpass filtering')
                            else
                                err.identifier
                                error('Have you added EP_DEN to path in Matlab?')
                            end
                        end
                    case 'NZT'
                        try
                            [coeff,denav,den_coeff,y,yo]= Run_NZT(av,stim,sc);
                        catch err         
                            if strcmp(err.identifier, 'MATLAB:TooManyInputs')
                                error('You probably ran the GUI of EP_DEN and now Matlab confuses the two Run_NZT, closing and restarting Matlab at least helps')
                            elseif strcmp(err.identifier, 'MATLAB:nologicalnan')                                
                                error('Why actually? happened with too low high cutoff frequency for bandpass filtering')
                            else
                                err.identifier
                                error('Have you added EP_DEN to path in Matlab?')
                            end
                        end
                end
                YDEN=st_den(x,den_coeff,sc,samples,handles);
        
            case 'load_den_coeff'        
                
                 [coeff,denav,den_coeff,y,yo]= Run_NZT(av,handles);
                 [filename, pathname] = uigetfile('*.mat','Select file');
                 matfile=load([pathname filename]);
                 den_coeff=matfile.den_coeff;
                 [denav,y,den_coeff]=st_den(av,den_coeff,handles);
                 [YDEN]=st_den(x,den_coeff,sc,samples,handles);  
        end

        %% Plotting        
        if i == 1
            % plot_denSingleTrial(samples, stim, sr, av, denav, plot_type, sc, coeff, den_coeff, yo, y, YDEN, sweeps, max_contour)
        end
        
        %% Return the denoised ERP
        YDEN = YDEN'; % 32x2048 -> 2048x32, transpose matrix
        epoch = YDEN(:);% (32*2048)x1, into one vector as it came in

       