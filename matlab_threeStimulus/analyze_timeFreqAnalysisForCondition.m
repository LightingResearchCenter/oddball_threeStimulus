function [realCoefs, imagCoefs, TFStat] = analyze_timeFreqAnalysisForCondition(epochsIn, artifactIndices, erpType, IAF_peak, parameters, handles)

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempTimeFreq.mat';
        if nargin == 0
            load('debugPath.mat')
            load(fullfile(path.debugMATs, debugMatFileName))
            close all
        else
            if handles.flags.saveDebugMATs == 1
                if 1 == 1 % ~strcmp(erpType, 'standard')
                    % do not save for standard tone as there are so many
                    % trials that debugging and developing of this function
                    % is so much slower compared to target and distracter
                    path = handles.path;
                    save('debugPath.mat', 'path')
                    save(fullfile(path.debugMATs, debugMatFileName))            
                end        
            end
        end 
    end

    numberOfEpochsIn = length(epochsIn.ERP);        
    [noOfSamples, noOfChannels] = size(epochsIn.ERP{1});   
    
    handles.parameters.timeFreq.Fc = 1; % central frequency?
        % restriction (δ, same as bandwidth parameter), delta, see Peng et al. (2012), pg. 8, point (1)
        % http://dx.doi.org/10.1371/journal.pone.0034163
        
    if parameters.timeFreq.logFrq == 1
        freqs = logspace(parameters.timeFreq.scaleLimits(1), parameters.timeFreq.scaleLimits(2), parameters.timeFreq.numberOfFreqs);
        scales = handles.parameters.timeFreq.Fc * parameters.EEG.srate ./ freqs;
    else
        freqs = linspace(parameters.timeFreq.scaleLimits(1), parameters.timeFreq.scaleLimits(2), parameters.timeFreq.numberOfFreqs);
        scales = handles.parameters.timeFreq.Fc * parameters.EEG.srate ./ freqs;
    end
    
    % original time vector in, e.g. -500 ms to 500 ms
    timeVectorIn = (linspace(-parameters.oddballTask.ERP_baseline, parameters.oddballTask.ERP_duration, noOfSamples))';

    % reduce the time resolution so we get some frequency resolution as
    % well divided by the downsampling factor "parameters.timeFreq.timeResolutionDivider"
    % fastwavelet.m requires the time points as indices            
    points = round(linspace(1, noOfSamples, noOfSamples/parameters.timeFreq.timeResolutionDivider));
        % round to make sure that there are no non-integers being
        % passed as input arguments that will be used as indices
        % inside the fastwavelet()    
        
    %% Compute Time-Freq           
    
        EEG_AllTheEpochs = zeros(size(epochsIn.ERP{1},1), parameters.EEG.nrOfChannels, numberOfEpochsIn);    
        disp(['      .. combining epochs to a matrix (', erpType, ')'])
        for epoch = 1 : numberOfEpochsIn            
            for ch = 1 : parameters.EEG.nrOfChannels                
                if artifactIndices(epoch, ch) == 1
                    EEG_AllTheEpochs(:,ch,epoch) = NaN;
                else
                    EEG_AllTheEpochs(:,ch,epoch) = epochsIn.ERP{epoch}(:,ch);
                end                
            end            
            % vectorize maybe this later                       
                % timep and freq should be the same for all epochs
                % disp(['CH: ', (num2str(ch))])
                % whos            
        end 
        
        [realCoefs, imagCoefs, realCoefs_SD, imagCoefs_SD, timep, freq, isNaN] = ...
                analyze_timeFreqWrapper(EEG_AllTheEpochs, parameters, epoch, scales, points, timeVectorIn, erpType, IAF_peak, handles);
        
        % Assign to output structure
        TFStat.real.mean = realCoefs;
        TFStat.real.SD = realCoefs_SD;
        TFStat.imag.mean = imagCoefs;
        TFStat.imag.SD = imagCoefs_SD;
        TFStat.n = sum(isNaN == 0)
            
    
    function [realStat, imagStat, n] = analyze_timeFreqComputeStats(realCoefs, imagCoefs, isNaN, dim, parameters, handles)
                
        [numberOfEpochsIn, noOfFreqs, noOfTimepoints, noOfChannels] = size(realCoefs);
        n = length(isNaN(isNaN == 0)); % number of non NaNs (i.e. non-artifacted epochs)
                
        realStat.mean = squeeze(nanmean(realCoefs, dim));
        realStat.SD = squeeze(nanstd(realCoefs, 0, dim));        
        imagStat.mean = squeeze(nanmean(imagCoefs, dim));
        imagStat.SD = squeeze(nanstd(imagCoefs, 0, dim));
        
            
        
        