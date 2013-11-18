function matrixStat = batch_calculateStatsForMatrix(matrixPerChannel, dim, flag, handles)

    [noOfTrials, noOfSubjects] = size(matrixPerChannel);
    matrix2D = matrixPerChannel;
    
    % matrixPerChannel
    % whos

    % remove Inf from input data (and see why this actually occurs?)
    infIndices = isinf(matrix2D);
    matrix2D(infIndices) = NaN;

    % average
    matrixStat.mean = nanmean(matrix2D, dim);
    matrixStat.SD = nanstd(matrix2D, flag, dim);
    matrixStat.n = (sum(~isnan(matrix2D')))';

    % a = matrixStat.mean
        % b = matrixStat.SD