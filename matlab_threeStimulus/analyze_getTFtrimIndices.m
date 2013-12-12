function indicesTrim = analyze_getTFtrimIndices(varToTrim, trimLimits, powerRaw, coiRaw, parameters)

    % low cut
    indexTemp = find(varToTrim == trimLimits(1));           
    if isempty(~indexTemp)      
        freqDiff = (varToTrim - trimLimits(1));
        freqDiff(freqDiff > 0) = NaN;
        [val1, fixed_ind1] = min(freqDiff);
        %disp(['    The closest match is = ', num2str(varToTrim(fixed_ind1)), ' Hz (index = ', num2str(fixed_ind1), ')'])
        indexTemp = fixed_ind1;
    end
    indicesTrim(1) = indexTemp; % only one value should be found

    % hi cut
    indexTemp = find(varToTrim == trimLimits(2));           
    if isempty(~indexTemp)      
        freqDiff = (varToTrim - trimLimits(2));
        % get the minimum of values above zero to ensure that the
        % frequency picked is higher than the high cutoff
        freqDiff(freqDiff < 0) = NaN;
        [val2, fixed_ind2] = min(freqDiff);
        %disp(['    The closest match is = ', num2str(varToTrim(fixed_ind2)), ' Hz (index = ', num2str(fixed_ind2), ')'])
        indexTemp = fixed_ind2;
    end
    indicesTrim(2) = indexTemp; % only one value should be found