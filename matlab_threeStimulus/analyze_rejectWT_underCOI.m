function [power, nanMask] = analyze_rejectWT_underCOI(time, freq, power, coi, parameters, handles)

    % vectorize later
    rows = length(freq);
    cols = length(time);

    invCoi = 1 ./ coi;

    index = zeros(length(time),1);
    try
        for t = 1 : cols          
            % find the largest frequency before COI for given time point
            index(t) = find(freq < invCoi(t), 1, 'last');
        end
    catch err
        err
        %invCoi(t)
        %freq
        find(freq < invCoi(t), 1, 'last')
        errorStr = sprintf('%s\n%s\n%s', ['No frequencies found to reject, this means probably that you use too high temporal resolution so that the frequency resolution is bad'], ...
                            ['default value for parameters.timeFreq.timeResolutionDivider is 16 which seem to give okay resolutin, if you go higher you get coarser time resolution'], ...
                            ['now you are using value = ', num2str(parameters.timeFreq.timeResolutionDivider), ' -> timeRes: ', num2str(time(2)-time(1)*1000), ' ms']);
        error(errorStr)
    end

    % create NaN-mask
    nanMask = false(rows,cols);
    for ind = 1 : cols
        nanMask(1:index(ind), ind) = 1;
    end               

    power(nanMask) = NaN;

