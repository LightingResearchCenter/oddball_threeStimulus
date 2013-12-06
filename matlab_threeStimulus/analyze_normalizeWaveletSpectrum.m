function powerDownsampled = analyze_normalizeWaveletSpectrum(powerDownsampled, t, f, baselineLimits, timeResolutionDivider, parameters, handles)

    %{
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempTFnorm.mat';
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

    timeRes = t(2) - t(1);
    freqRes = f(2) - f(1);   
    sampleRate = parameters.EEG.srate / timeResolutionDivider;
    %}  
    
    %baselineLimits

    % find indices corresponding to the baseline time
    indicesTrim = analyze_getTFtrimIndices(t, baselineLimits, [], [], parameters);
    
    % get the baseline mean (one mean for each frequency)
    baselineVector = nanmean(powerDownsampled(:,indicesTrim(1):indicesTrim(2)),2);
    
    % replicate the vector to a matrix so that the matrix dimensions match
    baselineVector = repmat(baselineVector, 1, size(powerDownsampled,2));
    
    % normalize so that we get a percentage difference like in the Figure 3
    % of Peng et al. (2012), http://dx.doi.org/10.1371/journal.pone.0034163
    % for example
    diff = powerDownsampled - baselineVector;
    powerDownsampled = 100 * (diff ./ (baselineVector));
    
    
    