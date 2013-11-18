function [peak, gravity] = analyze_findIndividualAlphaPeak(f, powerSpectrum, alphaRange)

    %% From Klimesh 1999

    % Inputs

    %    f              - frequency vector
    %    powerSpectrum  - corresponding power spectrum (of f)
    %    alphaRange     - your "ROI", e.g. [7 13] 
    
    % Klimesch W. 1999. 
    % EEG alpha and theta oscillations reflect cognitive and memory performance: a review and analysis. 
    % Brain Research Reviews 29:169–195. 
    % http://dx.doi.org/10.1016/S0165-0173(98)00056-3.
        
        % Usually, alpha frequency is defined in terms of peak or gravity frequency within the traditional 
        % alpha frequency range (f1 to f2) of about 7.5–12.5 Hz. Peak frequency is that spectral component 
        % within f1 to f2 which shows the largest power estimate (cf. Fig. 1A). Alpha frequency can also be 
        % calculated in terms of gravity (or `mean') frequency which is the weighted sum of spectral estimates, 
        % divided by alpha power: (∑(a(f)×f))/(∑ a(f)). Power spectral estimates at frequency f are denoted a(f). 
        % The index of summation is in the range of f1 to f2. Particularly if there are multiple peaks in the 
        % alpha range (for a classification see e.g., Ref. [39]), 
        
        % Klimesch defines the alpha band to have a range of (IAF - 3.5 to
        % 4 Hz --> IAF --> IAF + 1-1.5 Hz), e.g. with an IAF of 11 Hz, the
        % alpha band would be 7 - 12.5 Hz with the highest spread)
        
        %% gravity frequency appears the more adequate estimate of alpha frequency.


    %% Find indices corresponding to the desired alpha range, e.g. [7 13]
        ind1 = find(f == alphaRange(1));
        ind2 = find(f == alphaRange(2));
        
            % add some fix later if the exact freq is not found from input
            % f vector

    %% get the peak            
        [peakValue, peakInd] = max(powerSpectrum(ind1:ind2));
        peakInd = peakInd + ind1; % add the offset of ind1            
        peak = f(peakInd);
        
            if isnan(peakValue)
                peak = NaN;
            end

    %% calculate the gravity        
    
        SpectralEstimates = powerSpectrum(ind1:ind2) .* f(ind1:ind2);
        sumOfSpectralEstimates = nansum(SpectralEstimates);            
        powerSum = nansum(powerSpectrum(ind1:ind2));

        % gravity now is the sum of spectral estimates divided by the power
        gravity = sumOfSpectralEstimates / powerSum;


    