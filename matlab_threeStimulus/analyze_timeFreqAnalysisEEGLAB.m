function [ersp,itc,powbase,times,freqs,erspboot,itcboot] = analyze_timeFreqAnalysisEEGLAB(epochs_concan, erpType, parameters, handles)

    debugMatFileName = 'tempTimeFreqEEGLAB.mat';
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

    dataIn = epochs_concan.ERP;
    [numberOfSamples, numberOfChannels] = size(dataIn)
    
    epochLength = abs(-parameters.oddballTask.ERP_baseline - parameters.oddballTask.ERP_duration) * parameters.EEG.srate
    
    frames = epochLength
    tlimits = 1000*[-parameters.oddballTask.ERP_baseline parameters.oddballTask.ERP_duration]
    srate = parameters.EEG.srate
    
    for ch = 1 : numberOfChannels
        [ersp,itc,powbase,times,freqs,erspboot,itcboot] = timef(dataIn(:,ch)',frames,tlimits,srate)
    end
    