function epochsOut = pre_convertEpochs_to_NanVectors(epochsIn, artifactIndices)

    debugMatFileName = 'tempFASTERexcluding.mat';
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

    epochsOut = epochsIn;