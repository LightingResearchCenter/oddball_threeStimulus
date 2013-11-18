function powers = analyze_getPowersOfPredefinedBands(f, amplit, PSD, alpha, powerAnalysis, handles)
    
    debugMatFileName = 'tempPowerPredefinedBands.mat';
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
    
    IAF.amplit = alpha.amplit_gravity;
    IAF.PSD = alpha.PSD_gravity;
    
    powers = cell(length(powerAnalysis.eegBins.freqs), 1);
    
    for i = 1 : length(powerAnalysis.eegBins.freqs)
        
        % re-assign input to output
        powers{i}.label = powerAnalysis.eegBins.label{i};
        powers{i}.ch = powerAnalysis.eegBins.ch{i};
        powers{i}.freqVector = f;        
        
        disp(['            -- LABEL: ', powers{i}.label])
        
        %% Amplitude Spectrum
        
        ps_mean = (amplit{powerAnalysis.eegBins.ch{i}(1)}.mean + amplit{powerAnalysis.eegBins.ch{i}(2)}.mean) / 2;
        ps_SD = (amplit{powerAnalysis.eegBins.ch{i}(1)}.SD + amplit{powerAnalysis.eegBins.ch{i}(2)}.SD) / 2;
        
        if strcmp(powers{i}.label, 'lowAlpha') || strcmp(powers{i}.label, 'highAlpha') || strcmp(powers{i}.label, 'alpha')
            IAF_range = powerAnalysis.eegBins.freqs{i};
            powers{i}.ps.freqs = powerAnalysis.eegBins.freqs{i} + IAF.amplit;
            freqRange = powerAnalysis.alphaRange;
            [~, powers{i}.ps.powerData] = analyze_getBandAmplitude(f, ps_mean, ps_SD, freqRange, 'Power', IAF.amplit, IAF_range, handles);            
        else
            IAF_range = [4 4]; % not really used, just a dummy to avoid crashing
            powers{i}.ps.freqs = powerAnalysis.eegBins.freqs{i};
            freqRange = powers{i}.ps.freqs;
            [powers{i}.ps.powerData, ~] = analyze_getBandAmplitude(f, ps_mean, ps_SD, freqRange, 'Power', IAF.amplit, IAF_range, handles);            
        end
        
        %% PSD
        
        ps_mean = (PSD{powerAnalysis.eegBins.ch{i}(1)}.mean + PSD{powerAnalysis.eegBins.ch{i}(2)}.mean) / 2;
        ps_SD = (PSD{powerAnalysis.eegBins.ch{i}(1)}.SD + PSD{powerAnalysis.eegBins.ch{i}(2)}.SD) / 2;
        
        if strcmp(powers{i}.label, 'lowAlpha') || strcmp(powers{i}.label, 'highAlpha') || strcmp(powers{i}.label, 'alpha')
            IAF_range = powerAnalysis.eegBins.freqs{i};
            powers{i}.freqs = powerAnalysis.eegBins.freqs{i} + IAF.PSD;
            freqRange = powerAnalysis.alphaRange;
            [~, powers{i}.PSD.powerData] = analyze_getBandAmplitude(f, ps_mean, ps_SD, freqRange, 'PSD', IAF.PSD, IAF_range, handles);
        else
            IAF_range = [4 4]; % not really used, just a dummy to avoid crashing
            powers{i}.freqs = powerAnalysis.eegBins.freqs{i};
            freqRange = powers{i}.freqs;
            [powers{i}.PSD.powerData, ~] = analyze_getBandAmplitude(f, ps_mean, ps_SD, freqRange, 'PSD', IAF.PSD, IAF_range, handles);
            
        end
        
        % debug
        % psd = powers{i}.PSD.powerData
        
    end
    
    % make sure that the computation makes sense
    
        % Quick'n'dirty, fix later maybe, does not work if you add more
        % bands in init_DefaultParameters as indices are hard-coded
        beta = powers{4}.PSD.powerData;
        theta = powers{5}.PSD.powerData;
        delta = powers{6}.PSD.powerData; 
        alpha = powers{1}.PSD.powerData;
        total_FzCz = powers{7}.PSD.powerData;
        total_OzPz = powers{8}.PSD.powerData;
        FzCz_denom = beta + theta + delta;
        PzOz_denom = alpha;
                
        if (FzCz_denom / total_FzCz) > 1            
            warning('Ratios of beta+theta+delta should not exceed the total power in any conditions so you have a bug in the code probably?')
        end
        
        if (PzOz_denom / total_OzPz) > 1            
            warning('Ratio of alpha should not exceed the total power in any conditions so you have a bug in the code probably?')
        end
    
    