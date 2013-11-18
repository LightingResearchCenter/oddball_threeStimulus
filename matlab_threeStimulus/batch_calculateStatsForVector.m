function statsOut = batch_calculateStatsForVector(dataVector, dim, handles)

    % should have as many values as there are trials, i.e. 40 for
    % target
    
    % dataVector
    % whos

    statsOut.mean = nanmean(dataVector, dim);
    statsOut.SD = nanstd(dataVector, flag, dim);
    statsOut.n = mean(sum(~isnan(dataVector)));