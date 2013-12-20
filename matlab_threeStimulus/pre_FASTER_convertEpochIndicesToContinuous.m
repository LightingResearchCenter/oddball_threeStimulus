function artifactIndicesOut = pre_FASTER_convertEpochIndicesToContinuous(artifactIndices, epochLength, parameters, handles)

    %{
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempConvertIndices.mat';
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
    end
    %}

    
    [noOfEpochs, noOfChannels] = size(artifactIndices);  
    artifactIndicesOut = false(noOfEpochs*epochLength, noOfChannels);
    
    for ep = 1 : noOfEpochs       
        for ch = 1 : noOfChannels
            
            ind1 = (ep - 1)*epochLength + 1;
            ind2 = ind1 + epochLength - 1;
            
            if artifactIndices(ep,ch) == 1
                artifactIndicesOut(ind1:ind2, ch) = 1;
            else
                % nothing
            end
            
        end        
    end
    
    % vectorize maybe later for speed