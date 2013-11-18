function [windowType, r, segmentLength, nfft, freqRange, nOverlap, dataOut] = debug_setPowerAnalysisDefaults(dataIn, Fs, parameters)
        
    segmentLength = parameters.powerAnalysis.segmentLength * Fs;        
    nfft = segmentLength;
    nOverlap = parameters.powerAnalysis.nOverlap; % in %        
    r = parameters.powerAnalysis.tukeyWindowR;       
    windowType = 'Tukey';
    freqRange = 0 : 0.1 : parameters.filter.bandPass_hiFreq;

    % trim the channel data to be an integer of sampling frequency        
    numberOfSegmentLengths = floor(length(dataIn) / segmentLength);
    dataOut = dataIn(1:(numberOfSegmentLengths*segmentLength));
