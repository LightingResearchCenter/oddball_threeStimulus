function statsOut = batch_calculateStatsForVector(dataVector, dim, handles)
    
    % no of data points, e.g. number of subjects
    % [noOfSessions, NoOfDataPoints] = size(dataVector)
    flag = 0;
    statsOut.median = nanmean(dataVector, dim);
    statsOut.mean = nanmean(dataVector, dim);
    statsOut.SD = nanstd(dataVector, flag, dim);
    statsOut.n = mean(sum(~isnan(dataVector)));