function [rrTimes, rrPeakInterval, rrPeakAmplitude] = pre_createRRintervalVector(rr, rs, parameters, handles)

    numberOfSamples = length(rr);
    
    rrPeakInterval = zeros(length(rr)-1,1);
    rrPeakAmplitude = zeros(length(rr)-1,1);
    rrTimes = zeros(length(rr)-1,1);
    
    for i = 2 : numberOfSamples
        rrPeakInterval(i-1) = rr(i) - rr(i-1);
        rrPeakAmplitude(i-1) = rs(i);        
        rrTimes(i-1) = rr(i);
    end