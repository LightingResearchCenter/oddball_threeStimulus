function epochs_corr = pre_correctEpochByEpoch(epochsIn, dataType, parameters, handles)

    debugMatFileName = 'tempEpochByEpoch.mat';
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

    % assign RTs directly to output
    epochs_corr.RT_regular = epochsIn.RT_regular;
    epochs_corr.RT_irregular = epochsIn.RT_irregular;
    
    % also the samples per epoch
    epochs_corr.samplesPerEpoch = epochsIn.samplesPerEpoch;

    % REGULAR
    for i = 1 : length(epochsIn.oddball_regular)
        epochs_corr.oddball_regular{i} = pre_epochwiseCorrection(i, epochsIn.oddball_regular{i}, dataType, 'regular', parameters, handles);
    end


    % IRREGULAR
    for i = 1 : length(epochsIn.oddball_irregular)
        epochs_corr.oddball_irregular{i} = pre_epochwiseCorrection(i, epochsIn.oddball_irregular{i}, dataType, 'irregular', parameters, handles);
    end
   
    
    function epochCorr = pre_epochwiseCorrection(i, epochIn, dataType, condition, parameters, handles)
        
        % epochCorr
        % rows - as many samples as there is per epoch (i.e. 4096)
        % cols - as many channels as there (4 EEG channels + 1 EOG + 1 ECG)
        [rows,cols] = size(epochIn);
        
        % we could use the 'remove baseline' function from EEGLAB
        % Note that if data consists of multiple discontinuous epochs, each epoch should be separately baseline-zero'd using 
        % >> data = rmbase(data,frames,basevector) data = rmbase(data,frames,basevector);
        % http://sccn.ucsd.edu/eeglab/allfunctions/rmbase.m        
        % http://cognitrn.psych.indiana.edu/busey/temp/eeglabtutorial4.301/allfunctions/runica.html
        epochCorr = epochIn;
        
        % CHANNEL-WISE
        %{
        for i = 1 : cols
            epochCorr(:,i) = rmbase(epochIn(:,i));
        end
        %}
        
        % global mean
        [epochCorr(:,1:4), epochMean] = rmbase(epochIn(:,1:4));
        
        % Plot the results
        if handles.flags.showDebugPlots == 1
            scrsz = handles.style.scrsz;
            close all
            fig = figure('Color', 'w');
                set(fig, 'Position', [0.35*scrsz(3) 0.08*scrsz(4) 0.55*scrsz(3) 0.75*scrsz(4)])
                plot_epochByEpochCorrection(i, fig, epochIn, epochCorr, epochMean, dataType, condition, parameters, handles)
        end