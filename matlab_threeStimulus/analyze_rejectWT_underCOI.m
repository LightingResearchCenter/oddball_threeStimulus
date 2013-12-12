function power = analyze_rejectWT_underCOI(time, freq, power, coi, handles)

    % vectorize later

    rows = length(freq);
    cols = length(time);

    invCoi = 1 ./ coi;
    plot(time,invCoi)
    ylim([min(freq) max(freq)])

    index = zeros(length(time),1);
    for t = 1 : cols          
        % find the largest frequency before COI for given time point
        index(t) = find(freq < invCoi(t), 1, 'last');
    end

    % create NaN-mask
    nanMask = false(rows,cols);
    for ind = 1 : cols
        nanMask(1:index(ind), ind) = 1;
    end               

    power(nanMask) = NaN;

