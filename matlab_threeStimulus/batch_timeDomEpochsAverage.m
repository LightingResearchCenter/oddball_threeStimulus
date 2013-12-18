function [epochsProcessed, debug] = batch_timeDomEpochsAverage(epochsIn, thisFileMarkedAsOutlier, parameters, handles)

    %% DEBUG
    %{
    debugMatFileName = 'tempTimeDomainEpochsAverage.mat';
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
    %}
    
    % e.g.
    %{
    epochsIn = 

                    samplesPerEpoch: 5734
                          meanDelay: 596.8250
                       meanDelayStd: 94.5781
                  meanDelayMilliSec: 145.7092
               meanDelayMilliSecStd: 23.0904
               meanStimulusDuration: 878.9500
            meanStimulusDurationStd: 1.7967
       meanStimulusDurationMilliSec: 214.5874
    meanStimulusDurationMilliSecStd: 0.4387
                                ERP: {1x40 cell}
                                 RT: [1x40 double]
                         alphaFixed: [1x40 double]
                           alphaIAF: [1x40 double]
    %}
    
    % output debug info
    debug.meanDelay = epochsIn.meanDelay;
    debug.meanDelayStd = epochsIn.meanDelayStd;
    debug.meanStimulusDuration = epochsIn.meanStimulusDuration;
    debug.meanStimulusDurationStd = epochsIn.meanStimulusDurationStd;
        
    % get data dimensions
    noOfEpochsIn = length(epochsIn.ERP);
    [noOfSamplesIn, noOfChannels] = size(epochsIn.ERP{1});
    
    if thisFileMarkedAsOutlier == 0
    
        % Preallocate matrix 
        % (do not save EOG and ECG, thus the EEG.nrOfChannels
        allTheEpochs = zeros(noOfSamplesIn, parameters.EEG.nrOfChannels, noOfEpochsIn);
        isNaN = false(noOfEpochsIn, parameters.EEG.nrOfChannels);

        % go through the epochs   
        for ep = 1 : noOfEpochsIn 
            allTheEpochs(:,:,ep) = epochsIn.ERP{ep}(:,1:parameters.EEG.nrOfChannels);
            isNaN(ep,:) = logical(sum(isnan(epochsIn.ERP{ep}(:,1:parameters.EEG.nrOfChannels))));
        end

        %whos

        % compute the stats
        dim = 3;
        flag = 0;
        epochsProcessed.mean = nanmean(allTheEpochs,dim);
        epochsProcessed.median = nanmedian(allTheEpochs,dim);
        epochsProcessed.SD = nanstd(allTheEpochs,flag,dim);
        epochsProcessed.n = sum(~isNaN);

    else
        
        % just return NaN vectors if this specific file is marked as an
        % outlier
        epochsProcessed.mean = zeros(noOfSamplesIn, parameters.EEG.nrOfChannels) * NaN;
        epochsProcessed.median = zeros(noOfSamplesIn, parameters.EEG.nrOfChannels) * NaN;
        epochsProcessed.SD = zeros(noOfSamplesIn, parameters.EEG.nrOfChannels) * NaN;
        epochsProcessed.n = zeros(1, parameters.EEG.nrOfChannels) * NaN;
        
    end