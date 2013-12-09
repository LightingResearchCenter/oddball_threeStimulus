function [dataOut, NaN_indices, numberOfNaNs] = pre_artifactFixedThreshold(dataIn, EOG, parameters, handles)
    
    %{
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempArtifactRemoval.mat';
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
    

    [dataSamples, dataChannels] = size(dataIn);
    
    % Use short variable names for detrending filtering
    loFreq = parameters.artifacts.fixedDetrendingLowCut;
    hiFreq = parameters.artifacts.fixedDetrendingHighCut;
    N = parameters.artifacts.fixedDetrendingOrder;
    dataIn_detrended = zeros(dataSamples, dataChannels);

    %% Use fixed threshold        
    
        if parameters.artifacts.applyFixedThrRemoval == 1     
            
            disp(['    - ARTIFACT STATISTICS'])

            %{
            for i = 1 : dataChannels
                dataIn_DCoff(:,i) = (removedc(dataIn(:,i)', round(parameters.EEG.srate/2)))'; % remove DC offset first, function from ERPLAB
                dataIn_detrended(:,i) = detrend(dataIn_DCoff(:,i),'linear'); % detrend
            end
            %}
            % The above not really sufficient, better maybe to apply a
            % band-pass filter with a high cutoff frequency            

            for j = 1 : dataChannels                       
                dataIn_detrended(:,j) = pre_bandbassFilter(dataIn(:,j), parameters.EEG.srate, [hiFreq, loFreq], N, N*10, handles);   
                
                % check if output is valid
                nanIndices = isnan(dataIn_detrended(:,j));
                if length(dataIn_detrended(nanIndices,j)) == length(dataIn_detrended(:,j))
                    warning('NaN vector returned from bandpass filter, you probably used too high order, try to reduce and what happens, auto-reduce the order by 2 now')                    
                    N = N - 2; % reduce order, and try again
                    dataIn_detrended(:,j) = pre_bandbassFilter(dataIn(:,j), parameters.EEG.srate, [hiFreq, loFreq], N, N*10, handles);
                    nanIndices = isnan(dataIn_detrended(:,j));
                    if length(dataIn_detrended(nanIndices,j)) == length(dataIn_detrended(:,j))
                       error('NaN vector still returned, check what is the problem!') 
                    end                    
                end
            end

            % define indices
            NaN_indices_fixed = abs(dataIn_detrended) > parameters.artifacts.fixedThr; 

            % status on command window of how many artifacts were found
            numberOfNaNs = length(NaN_indices_fixed(NaN_indices_fixed == 1));
            NaNPercentage = (numberOfNaNs / (dataSamples*dataChannels)) * 100;
            disp(['       FIXED'])
            disp(['       .. Fixed threshold - ', num2str(parameters.artifacts.fixedThr), ' uV, ', 'Number of NaNs: ', num2str(numberOfNaNs), ', percentage: ', num2str(NaNPercentage), '%'])

        else
            NaN_indices_fixed = zeros(length(dataIn(:,1)),1);
            disp(['     .. Artifacts not searched with the fixed EEG threshold'])
        end
    
    %% Find excessive eye movements 
    
        if parameters.artifacts.applyFixedThrEOGRemoval == 1
    
            % Detrend EOG first        
            EOG = pre_bandbassFilter(EOG, parameters.EEG.srate, [hiFreq, loFreq], N, N*10, handles);                

            % Find indices of threshold exceeding voltages
            NaN_indices_EOG = abs(EOG) > parameters.artifacts.fixedThrEOG;

            % status on command window of how many artifacts were found
            numberOfNaNs_EOG = length(NaN_indices_EOG(NaN_indices_EOG == 1));
            NaNPercentage = (numberOfNaNs_EOG / dataSamples) * 100;
            disp(['       .. .. EOG threshold - ', num2str(parameters.artifacts.fixedThrEOG), ' uV, ',  'Number of NaNs: ', num2str(numberOfNaNs_EOG), ', percentage: ', num2str(NaNPercentage), '%'])
            
        else
            
            NaN_indices_EOG = zeros(length(dataIn(:,1)),1);
            disp(['     .. Artifacts not searched with the fixed EOG threshold'])  
            
        end
        
    %% C.R.A.P from ERPLAB
    
        
        % Use EEGLAB structure fields
        EEG.data = dataIn';
        EEG.epoch = [];
        EEG.srate = parameters.EEG.srate;
        EEG.pnts = length(dataIn);
        EEG.event(1).type = 0;
        EEG.xmin = 0;
        EEG.xmax = length(EEG.data) / parameters.EEG.srate;
        EEG.setname = 'dummy blank string';
        
        if parameters.artifacts.applyContinuousCRAP == 1            
            crap_NaNIndices = crap_mod(EEG, parameters.artifacts.CRAP.continuous_ampth, parameters.artifacts.CRAP.continuous_windowWidth, parameters.artifacts.CRAP.continuous_windowStep, ...
                        1:parameters.EEG.nrOfChannels);
            noOfNans = sum(crap_NaNIndices);
            NaNPercentage = (noOfNans / dataSamples) * 100;
            disp(['      C.R.A.P - ', num2str(parameters.artifacts.CRAP.continuous_ampth(1)), ' - ', num2str(parameters.artifacts.CRAP.continuous_ampth(2)), ' uV, '...
                         'wWidth = ', num2str(parameters.artifacts.CRAP.continuous_windowWidth), ', wStep = ', num2str(parameters.artifacts.CRAP.continuous_windowStep), ' ms, '...
                         'Number of NaNs: ', num2str(noOfNans), ', percentage: ', num2str(NaNPercentage), '%'])                               
        end

    %% Combine the indices
    
        % Now the _fixed indices have 4 columns (or the number of different
        % EEG channels that you have), whereas the EOG indices vector only
        % have one column so we need to replicate the EOG to match the
        % channels. We could do the vice versa as well
        [chRows, chCols] = size(NaN_indices_fixed);
        [eogRows, eogCols] = size(NaN_indices_EOG);
        repCount = chCols / eogCols;
        
        % replicate the rows
        NaN_indices_EOG = repmat(NaN_indices_EOG, 1, 4);
        
        NaN_indices = logical(NaN_indices_fixed + NaN_indices_EOG);
        dataOut = dataIn_detrended;
        
        %{
        doNotRemoveArtifacs = 1;
        if doNotRemoveArtifacs == 1
            disp(['       .. .. Artifacts were not converted to NaNs at this point'])
        else
            dataOut(NaN_indices) = NaN;
            disp(['       .. .. Artifacts were converted to NaNs'])
        end
        %}
        dataOut(NaN_indices) = NaN;
        
        
    % DEBUG
    %{
    subplot(3,1,1)
        plot(dataIn)
    subplot(3,1,2)
        plot(dataIn_detrended)
        ylim([-200 200])
    subplot(3,1,3)
        plot(dataOut)
        ylim([-200 200])
        
    drawnow
    pause
    %}