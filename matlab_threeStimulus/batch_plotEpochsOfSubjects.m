function batch_plotEpochsOfSubjects(epochsMatrix, erpType, erpFilterType, epochType, fileNameFields, statParam, handles)

    %% DEBUG
    debugMatFileName = 'tempBatchEpochPlots.mat';
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
    
    whos
    
    % get fieldnames IN
    conditions = fieldnames(epochsMatrix)
    sessions = fieldnames(epochsMatrix.(conditions{1}))
    erpResponses = fieldnames(epochsMatrix.(conditions{1}).(sessions{1}))
    erpTypes = fieldnames(epochsMatrix.(conditions{1}).(sessions{1}).(erpResponses{1}))
    filterTypes = fieldnames(epochsMatrix.(conditions{1}).(sessions{1}).(erpResponses{1}).(erpTypes{1}))
    statFields = fieldnames(epochsMatrix.(conditions{1}).(sessions{1}).(erpResponses{1}).(erpTypes{1}).(filterTypes{1}))
    
    exIn = epochsMatrix.(conditions{1}).(sessions{1}).(erpResponses{1}).(erpTypes{1}).(filterTypes{1})
    
    [noOfDataSamplesPerERP, noOfChannels, noOfSubjects] = size(epochsMatrix.(conditions{1}).(sessions{1}).(erpResponses{1}).(erpTypes{1}).(filterTypes{1}).mean)
    
    % Customize the color scale, use the distinguishable_colors from Matlab
    % FileExchange by Tim Holy, http://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
    % See also: http://blogs.mathworks.com/pick/2008/08/15/colors-for-your-multi-line-plots/
    figure
    
        ColorSet = distinguishable_colors(noOfSubjects);
        set(gca, 'ColorOrder', ColorSet);

        for ch = 1 : noOfChannels
            subplot(4,1,ch)
            plot(squeeze(exIn.(statParam)(:,ch,:)))
            hold on
        end

    figure
    
        ColorSet = distinguishable_colors(noOfSubjects);
        set(gca, 'ColorOrder', ColorSet);

        for ch = 1 : noOfChannels
            subplot(4,1,ch)
            plot(nanmean(exIn.(statParam)(:,ch,:),3))
            hold on
        end
