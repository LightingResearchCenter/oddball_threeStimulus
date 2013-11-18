function [realCoefs, imagCoefs, TFStat] = analyze_timeFreqAnalysisForCondition(epochsIn, erpType, parameters, handles)

    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempTimeFreq.mat';
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
        
    % Now there are 6 channels coming in as the EOG and ECGs are
    % not trimmed away so we only go through the EEG channels here  
    realCoefs = zeros(numberOfEpochsIn, length(scales), length(points), parameters.EEG.nrOfChannels);
    imagCoefs = zeros(numberOfEpochsIn, length(scales), length(points), parameters.EEG.nrOfChannels);
    isNaN = zeros(numberOfEpochsIn, 1);

        
    %% Compute Time-Freq   
        
        disp(['      .. combining epochs to a matrix (', erpType, ')'])
        for epoch = 1 : numberOfEpochsIn            
            EEG_AllTheEpochs(:,:,epoch) = epochsIn.ERP{epoch}(:,1:parameters.EEG.nrOfChannels);        
                       
            % timep and freq should be the same for all epochs
            % disp(['CH: ', (num2str(ch))])
            % whos            
        end 
        
        [realCoefs, imagCoefs, timep, freq, isNaN] = ...
                analyze_timeFreqWrapper(EEG_AllTheEpochs, parameters, epoch, scales, points, timeVectorIn, erpType, handles);        

        
    %% Average the epochs 
    
        disp(['        .. averaging epochs (', erpType, ')'])
        % now the input matrix is 4D and we want to average the epochs,
        % which means the dimension should be 1
        dim = 1;
        [TFStat.real, TFStat.imag, TFStat.n] = analyze_timeFreqComputeStats(realCoefs, imagCoefs, isNaN, dim, parameters, handles);        
        TFStat.timep = timep;
        TFStat.freq = freqs;
        TFStat.scale = scales;
        try
            TFStat.Fc = Fc;
        catch
            warning('fix Fc')
            TFStat.Fc = 1;
        end
        TFStat.bandwidth = parameters.timeFreq.bandwidthParameter;

        
    %% DEBUG Plot
    
        %{
        close all
        figure    
        contourLevels = 64;        
        subplot(2,1,1)
        contourf(TFStat.real.mean(:,:,1), 64, 'EdgeColor', 'none'); title('MEAN'); xlabel('Time'); ylabel('Hz');
        colorbar
        subplot(2,1,2)
        contourf(TFStat.real.SD(:,:,1), 64, 'EdgeColor', 'none'); title('SD'); xlabel('Time'); ylabel('Hz');
        colorbar
        %whos        
        %}
                
    
    function [realStat, imagStat, n] = analyze_timeFreqComputeStats(realCoefs, imagCoefs, isNaN, dim, parameters, handles)
                
        [numberOfEpochsIn, noOfFreqs, noOfTimepoints, noOfChannels] = size(realCoefs);
        n = length(isNaN(isNaN == 0)); % number of non NaNs (i.e. non-artifacted epochs)
                
        realStat.mean = squeeze(nanmean(realCoefs, dim));
        realStat.SD = squeeze(nanstd(realCoefs, 0, dim));        
        imagStat.mean = squeeze(nanmean(imagCoefs, dim));
        imagStat.SD = squeeze(nanstd(imagCoefs, 0, dim));
        
            